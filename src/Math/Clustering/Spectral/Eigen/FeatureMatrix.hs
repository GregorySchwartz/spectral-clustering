{- Math.Clustering.Spectral.Eigen.FeatureMatrix
Gregory W. Schwartz

Collects the functions pertaining to sparse spectral clustering.
-}

{-# LANGUAGE BangPatterns #-}

module Math.Clustering.Spectral.Eigen.FeatureMatrix
    ( B (..)
    , B1 (..)
    , B2 (..)
    , LabelVector (..)
    , spectral
    , spectralCluster
    , spectralClusterK
    , getB
    , b1ToB2
    , getSimilarityFromB2
    ) where

-- Remote
import Data.Bool (bool)
import Data.Function (on)
import Data.List (sortBy, maximumBy, transpose)
import Data.Maybe (fromMaybe)
import Safe (headMay)
import qualified AI.Clustering.KMeans as K
import qualified Data.Eigen.SparseMatrix as S
import qualified Data.Map.Strict as Map
import qualified Data.Vector as V
import qualified Data.Vector.Storable as VS
import qualified Data.Vector.Unboxed as U
import qualified Numeric.LinearAlgebra as H
import qualified Numeric.LinearAlgebra.Devel as H
import qualified Numeric.LinearAlgebra.SVD.SVDLIBC as SVD

-- Local

-- | Output vector containing cluster assignment (0 or 1).
type LabelVector = S.SparseMatrixXd
-- | B1 observation by feature matrix.
newtype B1 = B1 { unB1 :: S.SparseMatrixXd } deriving (Show)
-- | B2 term frequency-inverse document frequency matrix of B1.
newtype B2 = B2 { unB2 :: S.SparseMatrixXd } deriving (Show)
-- | Diagonal matrix from \(diag(B(B^{T}1))\).
newtype D  = D { unD :: S.SparseMatrixXd } deriving (Show)
-- | Matrix from \(D^{-1/2}B}\).
newtype C  = C { unC :: S.SparseMatrixXd } deriving (Show)
-- | Normed rows of B2. For a complete explanation, see Shu et al., "Efficient
-- Spectral Neighborhood Blocking for Entity Resolution", 2011.
newtype B  = B { unB :: S.SparseMatrixXd } deriving (Show)

-- | Assign close to 0 as 0.
epsilonZero :: Double -> Double
epsilonZero x = if abs x < 1e-12 then 0 else x

-- | Normalize the input matrix by column. Here, columns are features.
b1ToB2 :: B1 -> B2
b1ToB2 (B1 b1) =
    B2
        . S._imap (\ i j x
                 -> (log (fromIntegral n / (fromMaybe 0 $ dVec VS.!? j))) * x
                  )
        $ b1
  where
    dVec :: VS.Vector Double
    dVec = maybe (error "Cannot get number of non-zeros.") (VS.map fromIntegral)
         . S.innerNNZs
         . S.uncompress
         $ b1
    n = S.rows b1

-- | Euclidean norm each row.
b2ToB :: B2 -> B
b2ToB (B2 b2) =
    B
        . S._imap (\ i j x
                  -> x / (fromMaybe (error "Norm is 0.") $ eVec VS.!? i)
                  )
        $ b2
  where
    eVec :: VS.Vector Double
    eVec = VS.fromList . fmap S.norm . S.getRows $ b2

-- | Get the signed diagonal transformed B matrix.
bToD :: B -> D
bToD (B b) = D . S.diagCol 0 $ (S._map abs b) * ((S._map abs $ S.transpose b) * S.ones n)
  where
    n = S.rows b

-- | Get the matrix C as input for SVD.
bdToC :: B -> D -> C
bdToC (B b) (D d) = C $ (S._map (\x -> x ** (- 1 / 2)) d) * b

-- | Obtain the second largest value singular vector (or Nth) and E on of a
-- sparse matrix.
secondLeft :: Int -> Int -> S.SparseMatrixXd -> S.SparseMatrixXd
secondLeft n e m = S.transpose
                 . S.fromDenseList
                 . fmap H.toList
                 . drop (n - 1)
                 . H.toRows
                 . (\(!x, _, _) -> x)
                 . SVD.sparseSvd (e + (n - 1))
                 . H.mkCSR
                 . fmap (\(!i, !j, !x) -> ((i, j), x))
                 . S.toList
                 $ m

-- | Get the normalized matrix B from an input matrix where the features are
-- columns and rows are observations. Optionally, do not normalize.
getB :: Bool -> S.SparseMatrixXd -> B
getB True = b2ToB . b1ToB2 . B1
getB False = b2ToB . B2

-- | Returns the second left singular vector (or Nth) of a sparse spectral
-- process. Assumes the columns are features and rows are observations. B is the
-- normalized matrix (from getB). See Shu et al., "Efficient Spectral
-- Neighborhood Blocking for Entity Resolution", 2011.
spectral :: Int -> Int -> B -> S.SparseMatrixXd
spectral n e b = S._map epsilonZero . secondLeft n e . unC . bdToC b . bToD $ b

-- | Returns a vector of cluster labels for two groups by finding the second
-- left singular vector of a special normalized matrix. Assumes the columns are
-- features and rows are observations. B is the normalized matrix (from getB).
-- See Shu et al., "Efficient Spectral Neighborhood Blocking for Entity
-- Resolution", 2011.
spectralCluster :: B -> LabelVector
spectralCluster (B b)
  | S.rows b < 1  = S.fromDenseList [[]]
  | S.rows b == 1 = S.fromDenseList [[0]]
  | otherwise     = S.fromDenseList
                  . (fmap . fmap) (bool 0 1 . (>= 0))
                  . S.toDenseList
                  . spectral 2 1
                  $ B b

-- | Returns a vector of cluster labels for two groups by finding the largest
-- singular vectors and on of a special normalized matrix and running kmeans.
-- Assumes the columns are features and rows are observations. B is the
-- normalized matrix (from getB). See Shu et al., "Efficient Spectral
-- Neighborhood Blocking for Entity Resolution", 2011.
spectralClusterK :: Int -> Int -> B -> LabelVector
spectralClusterK e k (B b)
  | S.rows b < 1  = S.fromDenseList [[]]
  | S.rows b == 1 = S.fromDenseList [[0]]
  | otherwise     = consensusKmeans 100
                  . V.fromList
                  . fmap U.fromList
                  . concatMap S.toDenseList
                  . S.getRows
                  . S.fromCols
                  . fmap normNormalize
                  . S.getCols
                  . spectral 2 e
                  $ B b

-- | Consensus kmeans.
consensusKmeans :: Int -> V.Vector (U.Vector Double) -> LabelVector
consensusKmeans x vs = S.fromDenseList
                     . fmap ((:[]) . fromIntegral . mostCommon)
                     . transpose
                     . fmap kmeansFunc
                     $ [1 .. fromIntegral x]
  where
    kmeansFunc run =
      (\xs -> if headMay xs == Just 1 then fmap (bool 0 1 . (== 0)) xs else xs)
        . U.toList
        . K.membership
        . K.kmeansBy 2 vs id
        $ K.defaultKMeansOpts
            { K.kmeansMethod = K.Forgy
            , K.kmeansClusters = False
            , K.kmeansSeed = U.fromList [run]
            }

-- | Get the most common element of a list.
mostCommon :: (Ord a) => [a] -> a
mostCommon [] = error "Cannot find most common element of empty list."
mostCommon [x] = x
mostCommon xs = fst
               . maximumBy (compare `on` snd)
               . Map.toAscList
               . Map.fromListWith (+)
               . flip zip [1,1..]
               $ xs

-- | Normalize by the norm of a vector.
normNormalize :: S.SparseMatrixXd -> S.SparseMatrixXd
normNormalize xs = S._map (/ norm) xs
  where
    norm = S.norm xs

-- | Get the cosine similarity between two rows using B2.
getSimilarityFromB2 :: B2 -> Int -> Int -> Double
getSimilarityFromB2 (B2 b2) i j =
    (((S.getRow i b2) * (S.transpose $ S.getRow j b2)) S.! (0, 0))
        / (S.norm (S.getRow i b2) * S.norm (S.getRow j b2))
