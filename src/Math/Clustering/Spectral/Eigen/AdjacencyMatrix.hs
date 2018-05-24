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
    ) where

-- Remote
import Control.Monad (replicateM)
import Data.Bool (bool)
import Data.Function (on)
import Data.KMeans (kmeansGen)
import Data.List (sortBy)
import Data.Maybe (fromMaybe)
import System.Random.MWC (createSystemRandom, uniform)
import qualified Data.Eigen.SparseMatrix as S
import qualified Data.Eigen.SparseMatrix.Utility as S
import qualified Numeric.LinearAlgebra as H
import qualified Numeric.LinearAlgebra.Devel as H
import qualified Numeric.LinearAlgebra.SVD.SVDLIBC as SVD
import qualified Statistics.Quantile as Stat
import Debug.Trace

-- Local

type LabelVector     = S.SparseMatrixXd
type AdjacencyMatrix = S.SparseMatrixXd

-- | Returns the eigenvector with the second smallest eigenvalue of the
-- symmetric normalized Laplacian L. Computes real symmetric part of L, so
-- ensure the input is real and symmetric. Diagonal should be 0s for adjacency
-- matrix. Clusters the eigenvector using kmeans into k groups.
spectralClusterKNorm :: Int -> AdjacencyMatrix -> LabelVector
spectralClusterKNorm k = kmeansVec k . spectralNorm

-- | Returns the eigenvector with the second smallest eigenvalue of the
-- symmetric normalized Laplacian L. Computes real symmetric part of L, so
-- ensure the input is real and symmetric. Diagonal should be 0s for adjacency
-- matrix. Clusters the eigenvector by sign.
spectralClusterNorm :: AdjacencyMatrix -> LabelVector
spectralClusterNorm = S._map (bool 0 1 . (>= 0)) . spectralNorm

-- | Returns the eigenvector with the second smallest eigenvalue of the
-- symmetric normalized Laplacian L. Computes real symmetric part of L, so
-- ensure the input is real and symmetric. Diagonal should be 0s for adjacency
-- matrix. Uses I + Lnorm instead of I - Lnorm to find second largest singular
-- value instead of second smallest for Lnorm.
spectralNorm :: AdjacencyMatrix -> S.SparseMatrixXd
spectralNorm mat = secondLeft lNorm
  where
    lNorm    = i + (S.transpose invRootD * (mat * invRootD))
    invRootD = S.diagCol 0
             . S._map (\x -> if x == 0 then x else (1 / (x ** 2)))
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
            . (:[])
            . fmap snd
            . sortBy (compare `on` fst)
            . concatMap (\(c, xs) -> fmap (\(i, _) -> (i, c)) xs)
            . zip [0..] -- To get cluster id.
            . kmeansGen ((:[]) . snd) k
            . zip [0..] -- To keep track of index.
            . concat
            . S.toDenseList

-- | Obtain the second largest value singular vector of a sparse matrix.
secondLeft :: S.SparseMatrixXd -> S.SparseMatrixXd
secondLeft m = S.fromDenseList
             . (:[])
             . H.toList
             . last
             . H.toRows
             . (\(!x, _, _) -> x)
             . SVD.sparseSvd 2
             . H.mkCSR
             . fmap (\(!i, !j, !x) -> ((i, j), x))
             . S.toList
             $ m

-- | Obtain the degree matrix.
getDegreeMatrix :: AdjacencyMatrix -> S.SparseMatrixXd
getDegreeMatrix = S.diagRow 0 . getDegreeVector

-- | Obtain the degree vector.
getDegreeVector :: AdjacencyMatrix -> S.SparseMatrixXd
getDegreeVector = S.getRowSums
