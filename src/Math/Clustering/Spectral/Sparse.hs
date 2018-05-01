{- Math.Clustering.Spectral.Sparse
Gregory W. Schwartz

Collects the functions pertaining to sparse spectral clustering.
-}

{-# LANGUAGE BangPatterns #-}

module Math.Clustering.Spectral.Sparse
    ( B (..)
    , B2 (..)
    , spectral
    , spectralCluster
    , getB
    , b1ToB2
    , getSimilarityFromB2
    ) where

-- Remote
import Data.Bool (bool)
import qualified Data.Sparse.Common as S
import qualified Numeric.LinearAlgebra as H
import qualified Numeric.LinearAlgebra.Devel as H
import qualified Numeric.LinearAlgebra.SVD.SVDLIBC as SVD

-- Local

type LabelVector = S.SpVector Double
newtype B1 = B1 { unB1 :: S.SpMatrix Double } deriving (Show)
newtype B2 = B2 { unB2 :: S.SpMatrix Double } deriving (Show)
newtype D  = D { unD :: S.SpMatrix Double } deriving (Show)
newtype C  = C { unC :: S.SpMatrix Double } deriving (Show)
newtype B  = B { unB :: S.SpMatrix Double } deriving (Show)

-- | Normalize the input matrix by column. Here, columns are features.
b1ToB2 :: B1 -> B2
b1ToB2 (B1 b1) =
    B2
        . S.fromListSM (n, m)
        . fmap (\ (!i, !j, !x)
               -> (i, j, (log (fromIntegral n / (S.lookupDenseSV j dVec))) * x)
               )
        . S.toListSM
        $ b1
  where
    dVec :: S.SpVector Double
    dVec = S.vr
         . fmap (sum . fmap (\x -> if x > 0 then 1 else 0))
         . S.toRowsL -- faster than toColsL.
         . S.transposeSM
         $ b1
    n = S.nrows b1
    m = S.ncols b1

-- | Euclidean norm each row.
b2ToB :: B2 -> B
b2ToB (B2 b2) =
    B
        . S.fromListSM (n, m)
        . fmap (\(!i, !j, !x) -> (i, j, x / (S.lookupDenseSV i eVec)))
        . S.toListSM
        $ b2
  where
    eVec :: S.SpVector Double
    eVec = S.vr . fmap S.norm2 . S.toRowsL $ b2
    n = S.nrows b2
    m = S.ncols b2

-- | Find the Euclidean norm of a vector.
norm2 :: S.SpVector Double -> Double
norm2 = sqrt . sum . fmap (** 2)

-- | Get the diagonal transformed B matrix.
bToD :: B -> D
bToD (B b) = D
           . S.diagonalSM
           . flip S.extractCol 0
           $ b
       S.#~# ((S.transposeSM b) S.#~# (S.fromColsL [S.onesSV n]))
  where
    n = S.nrows b

-- | Get the matrix C as input for SVD.
bdToC :: B -> D -> C
bdToC (B b) (D d) = C $ (fmap (\x -> x ** (- 1 / 2)) d) S.#~# b

-- | Obtain the second left singular vector of a sparse matrix.
secondLeft :: S.SpMatrix Double -> S.SpVector Double
secondLeft m = S.sparsifySV
             . S.fromListDenseSV (S.nrows m)
             . H.toList
             . (!! 1)
             . H.toRows
             . (\(!x, _, _) -> x)
             . SVD.sparseSvd 2
             . H.mkCSR
             . fmap (\(!i, !j, !x) -> ((i, j), x))
             . S.toListSM
             $ m

-- | Get the normalized matrix B from an input matrix where the features are
-- columns and rows are observations. Optionally, do not normalize.
getB :: Bool -> S.SpMatrix Double -> B
getB True = b2ToB . b1ToB2 . B1
getB False = b2ToB . B2

-- | Returns the second left singular vector of a sparse spectral process.
-- Assumes the columns are features and rows are observations. B is the
-- normalized matrix (from getB). See Shu et al., "Efficient Spectral
-- Neighborhood Blocking for Entity Resolution", 2011.
spectral :: B -> S.SpVector Double
spectral b = secondLeft . unC $ c
  where
    d = bToD b
    c = bdToC b d

-- | Returns a vector of cluster labels for two groups by finding the second
-- left singular vector of a special normalized matrix. Assumes the columns are
-- features and rows are observations. B is the normalized matrix (from getB).
-- See Shu et al., "Efficient Spectral Neighborhood Blocking for Entity
-- Resolution", 2011.
spectralCluster :: B -> LabelVector
spectralCluster = S.sparsifySV . fmap (bool 0 1 . (>= 0)) . spectral

-- | Get the cosine similarity between two rows using B2.
getSimilarityFromB2 :: B2 -> Int -> Int -> Double
getSimilarityFromB2 (B2 b2) i j =
    S.dot (S.extractRow b2 i) (S.extractRow b2 j)
        / (S.norm2 (S.extractRow b2 i) * S.norm2 (S.extractRow b2 j))
