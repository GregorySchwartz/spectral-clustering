{- Math.Clustering.Spectral.Eigen.AdjacencyMatrix
Gregory W. Schwartz

Collects the functions pertaining to spectral clustering.
-}

{-# LANGUAGE BangPatterns #-}

module Math.Clustering.Spectral.Eigen.AdjacencyMatrix
    ( spectralClusterKNorm
    , spectralClusterNorm
    , spectralNorm
    , getDegreeMatrix
    , secondLeft
    , AdjacencyMatrix (..)
    , LabelVector (..)
    ) where

-- Remote
import Control.Monad (replicateM)
import Data.Bool (bool)
import Data.Function (on)
import Data.List (sortBy)
import Data.Maybe (fromMaybe)
import Safe (headMay)
import System.Random.MWC (createSystemRandom, uniform)
import qualified AI.Clustering.KMeans as K
import qualified Data.Eigen.SparseMatrix as S
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U
import qualified Numeric.LinearAlgebra as H
import qualified Numeric.LinearAlgebra.Devel as H
import qualified Numeric.LinearAlgebra.SVD.SVDLIBC as SVD
import qualified Statistics.Quantile as Stat

-- Local

type LabelVector     = S.SparseMatrixXd
type AdjacencyMatrix = S.SparseMatrixXd

-- | Returns the clustering of the eigenvectors with the second smallest
-- eigenvalues of the symmetric normalized Laplacian L. Computes real symmetric
-- part of L, so ensure the input is real and symmetric. Diagonal should be 0s
-- for adjacency matrix. Clusters the eigenvector using kmeans into k groups.
spectralClusterKNorm :: Int -> Int -> AdjacencyMatrix -> LabelVector
spectralClusterKNorm e k mat
  | S.rows mat < 1  = S.fromDenseList [[]]
  | S.rows mat == 1 = S.fromDenseList [[0]]
  | otherwise       = kmeansVec k . spectralNorm 1 e $ mat

-- | Returns the clustering of the eigenvectors with the second smallest
-- eigenvalues of the symmetric normalized Laplacian L. Computes real symmetric
-- part of L, so ensure the input is real and symmetric. Diagonal should be 0s
-- for adjacency matrix. Clusters the eigenvector by sign.
spectralClusterNorm :: AdjacencyMatrix -> LabelVector
spectralClusterNorm mat
  | S.rows mat < 1  = S.fromDenseList [[]]
  | S.rows mat == 1 = S.fromDenseList [[0]]
  | otherwise       = S.fromDenseList
                    . (fmap . fmap) (bool 0 1 . (>= 0))
                    . S.toDenseList
                    . spectralNorm 2 1
                    $ mat

-- | Returns the eigenvector with the second smallest eigenvalue (or N start)
-- and E on of the symmetric normalized Laplacian L. Computes real symmetric
-- part of L, so ensure the input is real and symmetric. Diagonal should be 0s
-- for adjacency matrix. Uses I + Lnorm instead of I - Lnorm to find second
-- largest singular value instead of second smallest for Lnorm.
spectralNorm :: Int -> Int -> AdjacencyMatrix -> S.SparseMatrixXd
spectralNorm n e mat
    | e < 1 = error "Less than 1 eigenvector chosen for clustering."
    | n < 1 = error "N < 1, cannot go before first eigenvector."
    | otherwise = secondLeft n e lNorm
  where
    lNorm    = i + (S.transpose invRootD * (mat * invRootD))
    invRootD = S.diagRow 0
             . S._map (\x -> if x == 0 then x else x ** (- 1 / 2))
             . getDegreeVector
             $ mat
    i        = S.ident . S.rows $ mat

-- | Second largest eigenvector. Unused, untested, don't use.
secondLargest :: S.SparseMatrixXd -> IO S.SparseMatrixXd
secondLargest a = do
    (first, firstVal) <- powerIt a

    let firstScale = S.scale firstVal (first * S.transpose first)
        b          = a - firstScale

    fmap fst $ powerIt b

-- | Rayleigh quotient. Takes a matrix and an eigenvector to find the
-- corresponding eigenvalue. Unused, untested, don't use.
rayQuot :: S.SparseMatrixXd -> S.SparseMatrixXd -> Double
rayQuot a x = ((S.transpose (a * x) * x) S.! (0, 0))
            / ((S.transpose x * x) S.! (0, 0))

-- | Power iteration. Unused, untested, don't use.
powerIt :: S.SparseMatrixXd -> IO (S.SparseMatrixXd, Double)
powerIt a = do
    g     <- createSystemRandom
    start <-
        fmap (S.fromDenseList . fmap (: [])) . replicateM (S.rows a) . uniform $ g
    let go
            :: Int
            -> Double
            -> S.SparseMatrixXd -- Eigenvector guess.
            -> Double -- Eigenvalue guess.
            -> (S.SparseMatrixXd, Double)
        go !i !e !b !lambda =
            if (abs (lambda' - lambda) < e) || (i < 0)
                then (b', lambda')
                else go (i - 1) e b' lambda'
          where
            absMat :: S.SparseMatrixXd -> S.SparseMatrixXd
            absMat = S._map abs
            b' :: S.SparseMatrixXd
            b' = S.scale (1 / maxElement ab) ab
            ab :: S.SparseMatrixXd
            ab = a * b
            maxElement = maximum . fmap (\(_, _, !x) -> abs x) . S.toList
            lambda' = S.norm ab / S.norm b
    return $ go 1000 0.000001 start (1 / 0)

-- | Executes kmeans to cluster a one dimensional vector.
kmeansVec :: Int -> S.SparseMatrixXd -> LabelVector
kmeansVec k = S.fromDenseList
            . fmap ((:[]) . fromIntegral)
            . U.toList
            . K.membership
            . (\x -> K.kmeansBy k x id K.defaultKMeansOpts)
            . V.fromList
            . fmap U.fromList
            . concatMap S.toDenseList
            . S.getRows
            . S.fromCols
            . fmap normNormalize
            . S.getCols

-- | Normalize by the norm of a vector.
normNormalize :: S.SparseMatrixXd -> S.SparseMatrixXd
normNormalize xs = S._map (/ norm) xs
  where
    norm = S.norm xs

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

-- | Obtain the second largest value singular vector of a sparse matrix.
denseSecondLeft :: S.SparseMatrixXd -> S.SparseMatrixXd
denseSecondLeft m = S.fromDenseList
                  . fmap (:[])
                  . H.toList
                  . (!! 2)
                  . H.toColumns
                  . (\(!x, _, _) -> x)
                  . H.svd
                  . H.assoc (S.rows m, S.cols m) 0
                  . fmap (\(!i, !j, !x) -> ((i, j), x))
                  . S.toList
                  $ m

-- | Obtain the signed degree matrix. Faster for columns.
getDegreeMatrix :: AdjacencyMatrix -> S.SparseMatrixXd
getDegreeMatrix = S.diagRow 0 . getDegreeVector

-- | Obtain the signed degree vector. Faster for columns.
getDegreeVector :: AdjacencyMatrix -> S.SparseMatrixXd
getDegreeVector = S.getColSums . S._map abs
