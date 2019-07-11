{- Math.Clustering.Spectral.Dense
Gregory W. Schwartz

Collects the functions pertaining to spectral clustering.
-}


{-# LANGUAGE BangPatterns #-}

module Math.Clustering.Spectral.Dense
    ( spectralClusterKNorm
    , spectralClusterNorm
    , spectralNorm
    , getDegreeMatrix
    , AdjacencyMatrix (..)
    , LabelVector (..)
    , B (..)
    , B1 (..)
    , B2 (..)
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
import qualified Data.Map.Strict as Map
import qualified Data.Vector as V
import qualified Data.Vector.Storable as VS
import qualified Data.Vector.Unboxed as U
import qualified Numeric.LinearAlgebra as H
import qualified Numeric.LinearAlgebra.Devel as H
import qualified Statistics.Quantile as S
import qualified Numeric.LinearAlgebra.SVD.SVDLIBC as SVD

-- Local

-- | Output vector containing cluster assignment (0 or 1).
type LabelVector     = H.Vector Double

-- | Adjacency matrix input.
type AdjacencyMatrix = H.Matrix Double

-- | B1 observation by feature matrix.
newtype B1 = B1 { unB1 :: H.Matrix Double } deriving (Show)
-- | B2 term frequency-inverse document frequency matrix of B1.
newtype B2 = B2 { unB2 :: H.Matrix Double } deriving (Show)
-- | Diagonal matrix from \(diag(B(B^{T}1))\). Here in a vector format.
newtype D  = D { unD :: H.Vector Double } deriving (Show)
-- | Matrix from \(D^{-1/2}B}\).
newtype C  = C { unC :: H.Matrix Double } deriving (Show)
-- | Normed rows of B2. For a complete explanation, see Shu et al., "Efficient
-- Spectral Neighborhood Blocking for Entity Resolution", 2011.
newtype B  = B { unB :: H.Matrix Double } deriving (Show)
-- | A diagonal matrix in a vector format.
newtype Diag  = Diag { unDiag :: H.Vector Double } deriving (Show)

-- | Assign close to 0 as 0.
epsilonZero :: Double -> Double
epsilonZero x = if abs x < 1e-12 then 0 else x

-- | Map hmatrix with indices.
cimap :: (Int -> Int -> Double -> Double) -> H.Matrix Double -> H.Matrix Double
cimap f mat = H.assoc (H.size mat) 0
            . concatMap (\ (!i, xs)
                        -> fmap (\ (!j, !x)
                                -> ( (i, j)
                                  , f i j x
                                  )
                                )
                                xs
                        )
            . zip [0..]
            . fmap (zip [0..])
            . H.toLists
            $ mat

-- | Normalize the input matrix by column. Here, columns are features.
b1ToB2 :: B1 -> B2
b1ToB2 (B1 b1) =
    B2
      . cimap (\ !i !j !x -> (log (fromIntegral n / (fromMaybe (error "Missing degree for observation. This would lead to divide by 0 error.") $ dVec VS.!? j))) * x)
      $ b1
  where
    dVec :: H.Vector Double
    dVec = H.fromList
         . fmap (H.sumElements . H.step)
         . H.toColumns
         $ b1
    n = H.rows b1
    m = H.cols b1

-- | Euclidean norm each row.
b2ToB :: B2 -> B
b2ToB (B2 b2) =
    B . cimap (\ !i !j !x -> x / (fromMaybe (error "Missing degree for observation. This would lead to divide by 0 error.") $ eVec VS.!? i)) $ b2
  where
    eVec :: H.Vector Double
    eVec = H.fromList . fmap H.norm_2 . H.toRows $ b2
    n = H.rows b2
    m = H.cols b2

-- | Get the signed diagonal transformed B matrix.
bToD :: B -> D
bToD (B b) = D
           . H.flatten
           $ (H.cmap abs b)
        H.<> ((H.cmap abs $ H.tr b) H.<> ((n H.>< 1) [1,1..]))
  where
    n = H.rows b

-- | Get the matrix C as input for SVD.
bdToC :: B -> D -> C
bdToC (B b) (D d) = C $ diagMatMult (Diag $ H.cmap (\x -> x ** (- 1 / 2)) d) b

-- | Obtain the second left singular vector (or N earlier) and E on of a sparse
-- matrix.
secondLeft :: Int -> Int -> H.Matrix Double -> [H.Vector Double]
secondLeft n e m =
  fmap (VS.drop (n - 1))
    . H.toColumns
    . (\(!x, _, _) -> x)
    . SVD.sparseSvd' (SVD.defaultSVDParams {SVD.maxIters = Just 1000}) (e + (n - 1))
    . H.mkCSR
    . filter (\((_, _), x) -> x /= 0)
    . concatMap (\(!i, xs) -> fmap (\(!j, !x) -> ((i, j), x)) xs)
    . zip [0..]
    . fmap (zip [0..])
    . H.toLists
    $ m

-- | Get the normalized matrix B from an input matrix where the features are
-- columns and rows are observations. Optionally, do not normalize.
getB :: Bool -> H.Matrix Double -> B
getB True = b2ToB . b1ToB2 . B1
getB False = b2ToB . B2

-- | Returns the second left singular vector (or from N) and E on of a sparse
-- spectral process. Assumes the columns are features and rows are observations.
-- B is the normalized matrix (from getB). See Shu et al., "Efficient Spectral
-- Neighborhood Blocking for Entity Resolution", 2011.
spectral :: Int -> Int -> B -> [H.Vector Double]
spectral n e b
    | e < 1     = error "Less than 1 eigenvector chosen for clustering."
    | n < 1 = error "N < 1, cannot go before first eigenvector."
    | otherwise =
        fmap (H.cmap epsilonZero) . secondLeft n e . unC . bdToC b . bToD $ b

-- | Returns a vector of cluster labels for two groups by finding the second
-- left singular vector of a special normalized matrix. Assumes the columns are
-- features and rows are observations. B is the normalized matrix (from getB).
-- See Shu et al., "Efficient Spectral Neighborhood Blocking for Entity
-- Resolution", 2011.
spectralCluster :: B -> LabelVector
spectralCluster (B b)
  | H.rows b < 1  = H.fromList []
  | H.rows b == 1 = H.fromList [0]
  | otherwise     = H.cmap (bool 0 1 . (>= 0))
                  . mconcat
                  . spectral 2 1
                  $ B b

-- | Returns a vector of cluster labels for two groups by finding the second
-- left singular vector and on of a special normalized matrix and running kmeans.
-- Assumes the columns are features and rows are observations. B is the
-- normalized matrix (from getB). See Shu et al., "Efficient Spectral
-- Neighborhood Blocking for Entity Resolution", 2011.
spectralClusterK :: Int -> Int -> B -> LabelVector
spectralClusterK e k (B b)
  | H.rows b < 1  = H.fromList []
  | H.rows b == 1 = H.fromList [0]
  | otherwise     = kmeansVec k . spectral 2 e $ B b

-- | Executes kmeans to cluster a vector.
kmeansVec :: Int -> [H.Vector Double] -> LabelVector
kmeansVec k = consensusKmeans 100
            . V.fromList
            . fmap U.convert
            . H.toRows
            . H.fromColumns
            . fmap H.normalize -- Normalize within eigenvectors (columns).
            . H.toColumns
            . H.fromRows
  where

-- | Consensus kmeans.
consensusKmeans :: Int -> V.Vector (U.Vector Double) -> LabelVector
consensusKmeans x vs = H.fromList
                     . fmap (fromIntegral . mostCommon)
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

-- | Get the cosine similarity between two rows using B2.
getSimilarityFromB2 :: B2 -> Int -> Int -> Double
getSimilarityFromB2 (B2 b2) i j =
    H.dot (H.flatten $ b2 H.? [i]) (H.flatten $ b2 H.? [j])
        / (H.norm_2 (H.flatten $ b2 H.? [i]) * H.norm_2 (H.flatten $ b2 H.? [j]))

-- | Returns the clustering of eigenvectors with the second smallest eigenvalues
-- and on of the symmetric normalized Laplacian L. Computes real symmetric part
-- of L, so ensure the input is real and symmetric. Diagonal should be 0s for
-- adjacency matrix. Clusters the eigenvector using kmeans into k groups from e
-- eigenvectors.
spectralClusterKNorm :: Int -> Int -> AdjacencyMatrix -> LabelVector
spectralClusterKNorm e k mat
  | H.rows mat < 1  = H.fromList []
  | H.rows mat == 1 = H.fromList [0]
  | otherwise       = kmeansVec k
                    . spectralNorm 2 e
                    $ mat

-- | Returns the eigenvector with the second smallest eigenvalue of the
-- symmetric normalized Laplacian L. Computes real symmetric part of L, so
-- ensure the input is real and symmetric. Diagonal should be 0s for adjacency
-- matrix.
spectralClusterNorm :: AdjacencyMatrix -> LabelVector
spectralClusterNorm mat
  | H.rows mat < 1  = H.fromList []
  | H.rows mat == 1 = H.fromList [0]
  | otherwise       = H.cmap (bool 0 1 . (>= 0))
                    . mconcat
                    . spectralNorm 2 1
                    $ mat

-- | Returns the eigenvectors with the Nth smallest eigenvalue and on of the
-- symmetric normalized Laplacian L. Computes real symmetric part of L, so
-- ensure the input is real and symmetric. Diagonal should be 0s for adjacency
-- matrix.
spectralNorm :: Int -> Int -> AdjacencyMatrix -> [H.Vector Double]
spectralNorm n e mat
    | e < 1 = error "Less than 1 eigenvector chosen for clustering."
    | n < 1 = error "N < 1, cannot go before first eigenvector."
    | otherwise = H.toRows
                . H.cmap epsilonZero -- Correct accuracy.
                . flip (H.??) (H.All, H.TakeLast e)
                . flip (H.??) (H.All, H.DropLast (n - 1))
                . snd
                . H.eigSH
                $ lNorm
  where
    lNorm = H.trustSym $ i - (matDiagMult (diagMatMult invD mat) invD)
    invD  = Diag
          . H.cmap (\x -> if x == 0 then x else x ** (- 1 / 2))
          . getDegreeVector
          $ mat
    i     = H.ident . H.rows $ mat

-- | Sort a matrix by values in a vector.
sortMatrixByVec :: H.Vector Double -> H.Matrix Double -> H.Matrix Double
sortMatrixByVec xs mat = mat H.¿ sortedIdx
  where
    sortedIdx =
      reverse . fmap fst . sortBy (compare `on` snd) . zip [0..] . H.toList $ xs

-- | Obtain the signed degree matrix.
getDegreeMatrix :: AdjacencyMatrix -> Diag
getDegreeMatrix = Diag . getDegreeVector

-- | Obtain the signed degree vector.
getDegreeVector :: AdjacencyMatrix -> H.Vector Double
getDegreeVector = H.vector . fmap (H.sumElements . H.cmap abs) . H.toRows

-- | Diagonal matrix multiply from vector type.
diagMatMult :: Diag -> H.Matrix Double -> H.Matrix Double
diagMatMult (Diag d) = cimap mult
  where
    mult i _ v = v * (d H.! i)

-- | Diagonal matrix multiply from vector type.
matDiagMult :: H.Matrix Double -> Diag -> H.Matrix Double
matDiagMult m (Diag d) = cimap mult m
  where
    mult _ j v = v * (d H.! j)
