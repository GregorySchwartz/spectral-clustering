{- Math.Clustering.Spectral.Dense
Gregory W. Schwartz

Collects the functions pertaining to spectral clustering.
-}

module Math.Clustering.Spectral.Dense
    ( spectralClusterKNorm
    , spectralClusterNorm
    , spectralNorm
    , spectralClusterP
    , spectralP
    , getDegreeMatrix
    ) where

-- Remote
import Data.Bool (bool)
import Data.Function (on)
import Data.KMeans (kmeansGen)
import Data.List (sortBy)
import qualified Numeric.LinearAlgebra as H
import qualified Statistics.Quantile as S

-- Local

type LabelVector     = H.Vector Double
type AdjacencyMatrix = H.Matrix Double

-- | Returns a vector of cluster labels by finding the eigenvector with the
-- largest eigenvalue of the random walk normalized Laplacian P. Computes real
-- symmetric part of P, so ensure the input is real and symmetric. Diagonal
-- should be 0s for adjacency matrix. Clusters the eigenvector using kmeans into
-- k groups.
spectralClusterKNorm :: Int -> AdjacencyMatrix -> LabelVector
spectralClusterKNorm k = H.fromList
                       . fmap snd
                       . sortBy (compare `on` fst)
                       . concatMap (\(c, xs) -> fmap (\(i, _) -> (i, c)) xs)
                       . zip [0..] -- To get cluster id.
                       . kmeansGen ((:[]) . snd) k
                       . zip [0..] -- To keep track of index.
                       . H.toList
                       . spectralNorm

-- | Returns a vector of cluster labels by finding the eigenvector with the
-- largest eigenvalue of the random walk normalized Laplacian P. Computes real
-- symmetric part of P, so ensure the input is real and symmetric. Diagonal
-- should be 0s for adjacency matrix.
spectralClusterNorm :: AdjacencyMatrix -> LabelVector
spectralClusterNorm = H.cmap (bool 0 1 . (>= 0)) . spectralNorm

-- | Returns a vector of cluster labels by finding the eigenvector with the
-- largest eigenvalue of the random walk normalized Laplacian P. Computes real
-- symmetric part of P, so ensure the input is real and symmetric. Diagonal
-- should be 0s for adjacency matrix.
spectralClusterP :: AdjacencyMatrix -> LabelVector
spectralClusterP mat = H.cmap (bool 0 1 . (> m)) eigVec
  where
    m = S.continuousBy S.s 2 4 eigVec
    eigVec = spectralP mat

-- | Returns the eigenvector with the largest eigenvalue of the random walk
-- normalized Laplacian P. Computes real symmetric part of P, so ensure the
-- input is real and symmetric. Diagonal should be 0s for adjacency matrix.
spectralP :: AdjacencyMatrix -> H.Vector Double
spectralP mat = head . H.toColumns . snd . H.eigSH $ rwLap
  where
    rwLap = H.sym $ invD H.<> mat
    invD  = H.diag
          . H.cmap (\x -> if x == 0 then x else (1 / x))
          . getDegreeVector
          $ mat

-- | Returns the eigenvector with the second smallest eigenvalue of the
-- symmetric normalized Laplacian L. Computes real symmetric part of L, so
-- ensure the input is real and symmetric. Diagonal should be 0s for adjacency
-- matrix.
spectralNorm :: AdjacencyMatrix -> H.Vector Double
spectralNorm mat = H.flatten
                 . flip (H.??) (H.All, H.PosCyc (H.idxs [-2]))
                 . snd
                 . H.eigSH
                 $ lNorm
  where
    lNorm = H.sym $ i - mconcat [invD, mat, invD]
    invD  = H.diag
          . H.cmap (\x -> if x == 0 then x else sqrt (1 / x))
          . getDegreeVector
          $ mat
    i     = H.ident . H.rows $ mat

-- | Obtain the degree matrix.
getDegreeMatrix :: AdjacencyMatrix -> H.Matrix Double
getDegreeMatrix = H.diag . getDegreeVector

-- | Obtain the degree vector.
getDegreeVector :: AdjacencyMatrix -> H.Vector Double
getDegreeVector = H.vector . fmap H.sumElements . H.toRows
