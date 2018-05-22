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
import Data.Bool (bool)
import Data.Function (on)
import Data.KMeans (kmeansGen)
import Data.List (sortBy)
import Data.Maybe (fromMaybe)
import Safe (headMay)
import qualified Data.Eigen.SparseMatrix as S
import qualified Data.Eigen.SparseMatrix.Utility as S
import qualified Numeric.LinearAlgebra as H
import qualified Numeric.LinearAlgebra.Devel as H
import qualified Numeric.LinearAlgebra.SVD.SVDLIBC as SVD
import qualified Statistics.Quantile as Stat

-- Local

type LabelVector     = S.SparseMatrixXd
type AdjacencyMatrix = S.SparseMatrixXd

-- | Returns a vector of cluster labels by finding the eigenvector with the
-- largest eigenvalue of the random walk normalized Laplacian P. Computes real
-- symmetric part of P, so ensure the input is real and symmetric. Diagonal
-- should be 0s for adjacency matrix. Clusters the eigenvector using kmeans into
-- k groups.
spectralClusterKNorm :: Int -> AdjacencyMatrix -> LabelVector
spectralClusterKNorm k = S.fromDenseList
                       . (:[])
                       . fmap snd
                       . sortBy (compare `on` fst)
                       . concatMap (\(c, xs) -> fmap (\(i, _) -> (i, c)) xs)
                       . zip [0..] -- To get cluster id.
                       . kmeansGen ((:[]) . snd) k
                       . zip [0..] -- To keep track of index.
                       . fromMaybe (error "More than one row in \"vector\".")
                       . headMay
                       . S.toDenseList
                       . spectralNorm

-- | Returns a vector of cluster labels by finding the eigenvector with the
-- largest eigenvalue of the random walk normalized Laplacian P. Computes real
-- symmetric part of P, so ensure the input is real and symmetric. Diagonal
-- should be 0s for adjacency matrix.
spectralClusterNorm :: AdjacencyMatrix -> LabelVector
spectralClusterNorm = S._map (bool 0 1 . (>= 0)) . spectralNorm

-- | Returns the eigenvector with the second smallest eigenvalue of the
-- symmetric normalized Laplacian L. Computes real symmetric part of L, so
-- ensure the input is real and symmetric. Diagonal should be 0s for adjacency
-- matrix.
spectralNorm :: AdjacencyMatrix -> S.SparseMatrixXd
spectralNorm mat = secondLeft lNorm
  where
    lNorm = i - (invD * mat * invD)
    invD  = S.diagCol 0
          . S._map (\x -> if x == 0 then x else sqrt (1 / x))
          . getDegreeVector
          $ mat
    i     = S.ident . S.rows $ mat

-- | Obtain the second lowest value singular vector of a sparse matrix
-- (different than with a feature matrix).
secondLeft :: S.SparseMatrixXd -> S.SparseMatrixXd
secondLeft m = S.fromDenseList
             . (:[])
             . H.toList
             . last
             . H.toRows
             . (\(!x, _, _) -> x)
             . SVD.sparseSvd (S.rows m - 1)
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
