{- Math.Clustering.Spectral.Dense
Gregory W. Schwartz

Collects the functions pertaining to spectral clustering.
-}

module Math.Clustering.Spectral.Dense
    ( spectralClusterKNorm
    , spectralClusterNorm
    , spectralNorm
    , getDegreeMatrix
    , AdjacencyMatrix (..)
    ) where

-- Remote
import Data.Bool (bool)
import Data.Function (on)
import Data.List (sortBy)
import Data.Maybe (fromMaybe)
import Safe (headMay)
import qualified AI.Clustering.KMeans as K
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U
import qualified Numeric.LinearAlgebra as H
import qualified Statistics.Quantile as S

-- Local

type LabelVector     = H.Vector Double
type AdjacencyMatrix = H.Matrix Double

-- | Returns the clustering of eigenvectors with the second smallest eigenvalues
-- and on of the symmetric normalized Laplacian L. Computes real symmetric part
-- of L, so ensure the input is real and symmetric. Diagonal should be 0s for
-- adjacency matrix. Clusters the eigenvector using kmeans into k groups from e
-- eigenvectors.
spectralClusterKNorm :: Int -> Int -> AdjacencyMatrix -> LabelVector
spectralClusterKNorm e k mat
  | H.rows mat < 1  = H.fromList []
  | H.rows mat == 1 = H.fromList [0]
  | otherwise       = V.convert
                    . U.map fromIntegral
                    . K.membership
                    . (\x -> K.kmeansBy k x id K.defaultKMeansOpts)
                    . V.fromList
                    . fmap V.convert
                    . H.toRows
                    . H.fromColumns
                    . fmap H.normalize -- Normalize within eigenvectors (columns).
                    . H.toColumns
                    . H.fromRows
                    . spectralNorm 1 e
                    $ mat

-- | Returns the eigenvector with the second smallest eigenvalue of the
-- symmetric normalized Laplacian L. Computes real symmetric part of L, so
-- ensure the input is real and symmetric. Diagonal should be 0s for adjacency
-- matrix.
spectralClusterNorm :: AdjacencyMatrix -> LabelVector
spectralClusterNorm mat
  | H.rows mat < 1  = H.fromList []
  | H.rows mat == 1 = H.fromList [0]
  | otherwise       =
      H.cmap (bool 0 1 . (>= 0)) . mconcat . spectralNorm 2 1 $ mat

-- | Returns the eigenvectors with the Nth smallest eigenvalue and on of the
-- symmetric normalized Laplacian L. Computes real symmetric part of L, so
-- ensure the input is real and symmetric. Diagonal should be 0s for adjacency
-- matrix.
spectralNorm :: Int -> Int -> AdjacencyMatrix -> [H.Vector Double]
spectralNorm n e mat
    | e < 1 = error "Less than 1 eigenvector chosen for clustering."
    | n < 1 = error "N < 1, cannot go before first eigenvector."
    | otherwise = H.toRows
                . flip (H.??) (H.All, H.TakeLast e)
                . flip (H.??) (H.All, H.DropLast (n - 1))
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
