{- SpectralClustering
Gregory W. Schwartz

Collects the functions pertaining to spectral clustering.
-}

module Math.SpectralClustering
    ( spectralClusterNorm
    , spectralNorm
    , spectralClusterP
    , spectralP
    , getDegreeMatrix
    ) where

-- Remote
import Data.Bool (bool)
import qualified Numeric.LinearAlgebra as H
import qualified Statistics.Quantile as S

-- Local

-- | Returns a vector of cluster labels by finding the eigenvector with the
-- largest eigenvalue of the random walk normalized Laplacian P. Computes real
-- symmetric part of P, so ensure the input is real and symmetric.
spectralClusterNorm :: H.Matrix Double -> H.Vector Double
spectralClusterNorm = H.cmap (bool 0 1 . (>= 0)) . spectralNorm

-- | Returns a vector of cluster labels by finding the eigenvector with the
-- largest eigenvalue of the random walk normalized Laplacian P. Computes real
-- symmetric part of P, so ensure the input is real and symmetric.
spectralClusterP :: H.Matrix Double -> H.Vector Double
spectralClusterP mat = H.cmap (bool 0 1 . (> m)) eigVec
  where
    m = S.continuousBy S.s 2 4 eigVec
    eigVec = spectralP mat

-- | Returns the eigenvector with the largest eigenvalue of the random walk
-- normalized Laplacian P. Computes real symmetric part of P, so ensure the
-- input is real and symmetric.
spectralP :: H.Matrix Double -> H.Vector Double
spectralP mat = head . H.toColumns . snd . H.eigSH $ rwLap
  where
    rwLap = H.sym $ H.inv d H.<> mat
    d     = getDegreeMatrix mat

-- | Returns the eigenvector with the second smallest eigenvalue of the
-- symmetric normalized Laplacian L. Computes real symmetric part of L, so
-- ensure the input is real and symmetric.
spectralNorm :: H.Matrix Double -> H.Vector Double
spectralNorm mat = head . drop 1 . reverse . H.toColumns . snd . H.eigSH $ lNorm
  where
    lNorm = H.sym $ i - mconcat [invD, mat, invD]
    invD  = H.inv . H.sqrtm $ d
    d     = getDegreeMatrix mat
    i     = H.ident . H.rows $ mat

-- | Obtain the degree matrix.
getDegreeMatrix :: H.Matrix Double -> H.Matrix Double
getDegreeMatrix = H.diag . H.vector . fmap H.sumElements . H.toRows
