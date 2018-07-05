{- Math.Clustering.Spectral.Dense
Gregory W. Schwartz

Collects the functions pertaining to spectral clustering.
-}

module Math.Clustering.Spectral.Dense
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
import qualified Numeric.LinearAlgebra as H
import qualified Statistics.Quantile as S

-- Local

type LabelVector     = H.Vector Double
type AdjacencyMatrix = H.Matrix Double

-- | Returns the eigenvector with the second smallest eigenvalue of the
-- symmetric normalized Laplacian L. Computes real symmetric part of L, so
-- ensure the input is real and symmetric. Diagonal should be 0s for adjacency
-- matrix. Clusters the eigenvector using kmeans into k groups.
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
-- | Returns the eigenvector with the second smallest eigenvalue of the
-- symmetric normalized Laplacian L. Computes real symmetric part of L, so
-- ensure the input is real and symmetric. Diagonal should be 0s for adjacency
-- matrix.
spectralClusterNorm :: AdjacencyMatrix -> LabelVector
spectralClusterNorm = H.cmap (bool 0 1 . (>= 0)) . spectralNorm

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
          . H.cmap (\x -> if x == 0 then x else x ** (- 1 / 2))
          . getDegreeVector
          $ mat
    i     = H.ident . H.rows $ mat

-- | Obtain the degree matrix.
getDegreeMatrix :: AdjacencyMatrix -> H.Matrix Double
getDegreeMatrix = H.diag . getDegreeVector

-- | Obtain the degree vector.
getDegreeVector :: AdjacencyMatrix -> H.Vector Double
getDegreeVector = H.vector . fmap H.sumElements . H.toRows
