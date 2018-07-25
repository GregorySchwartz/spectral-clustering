{- Math.Clustering.Spectral.Sparse
Gregory W. Schwartz

Collects the functions pertaining to sparse spectral clustering.
-}

{-# LANGUAGE BangPatterns #-}

module Math.Clustering.Spectral.Sparse
    ( B (..)
    , B1 (..)
    , B2 (..)
    , AdjacencyMatrix (..)
    , spectral
    , spectralCluster
    , spectralClusterK
    , spectralNorm
    , spectralClusterNorm
    , spectralClusterKNorm
    , getB
    , b1ToB2
    , getSimilarityFromB2
    ) where

-- Remote
import Data.Bool (bool)
import Data.Maybe (fromMaybe)
import Data.Function (on)
import Data.KMeans (kmeansGen)
import Data.List (sortBy, foldl1')
import Safe (headMay)
import qualified Data.Sparse.Common as S
import qualified Numeric.LinearAlgebra as H
import qualified Numeric.LinearAlgebra.Devel as H
import qualified Numeric.LinearAlgebra.SVD.SVDLIBC as SVD

-- Local

type LabelVector = S.SpVector Double
type AdjacencyMatrix = S.SpMatrix Double
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

-- | Obtain the second left singular vector and on of a sparse matrix.
secondLeft :: Int -> S.SpMatrix Double -> [S.SpVector Double]
secondLeft e m = fmap (S.sparsifySV . S.fromListDenseSV e . drop 1 . H.toList)
               . H.toColumns
               . (\(!x, _, _) -> x)
               . SVD.sparseSvd (e + 1)
               . H.mkCSR
               . fmap (\(!i, !j, !x) -> ((i, j), x))
               . S.toListSM
               $ m

-- | Get the normalized matrix B from an input matrix where the features are
-- columns and rows are observations. Optionally, do not normalize.
getB :: Bool -> S.SpMatrix Double -> B
getB True = b2ToB . b1ToB2 . B1
getB False = b2ToB . B2

-- | Returns the second left singular vector and on of a sparse spectral process.
-- Assumes the columns are features and rows are observations. B is the
-- normalized matrix (from getB). See Shu et al., "Efficient Spectral
-- Neighborhood Blocking for Entity Resolution", 2011.
spectral :: Int -> B -> [S.SpVector Double]
spectral e b
    | e < 1     = error "Less than 1 eigenvector chosen for clustering."
    | otherwise = secondLeft e . unC . bdToC b . bToD $ b

-- | Returns a vector of cluster labels for two groups by finding the second
-- left singular vector of a special normalized matrix. Assumes the columns are
-- features and rows are observations. B is the normalized matrix (from getB).
-- See Shu et al., "Efficient Spectral Neighborhood Blocking for Entity
-- Resolution", 2011.
spectralCluster :: B -> LabelVector
spectralCluster = S.sparsifySV
                . fmap (bool 0 1 . (>= 0))
                . foldl1' S.concatSV
                . spectral 1

-- | Returns a vector of cluster labels for two groups by finding the second
-- left singular vector and on of a special normalized matrix and running kmeans.
-- Assumes the columns are features and rows are observations. B is the
-- normalized matrix (from getB). See Shu et al., "Efficient Spectral
-- Neighborhood Blocking for Entity Resolution", 2011.
spectralClusterK :: Int -> Int -> B -> LabelVector
spectralClusterK e k = kmeansVec k . spectral e

-- | Executes kmeans to cluster a vector.
kmeansVec :: Int -> [S.SpVector Double] -> LabelVector
kmeansVec k = S.sparsifySV
            . S.vr
            . fmap snd
            . sortBy (compare `on` fst)
            . concatMap (\(c, xs) -> fmap (\(i, _) -> (i, c)) xs)
            . zip [0..] -- To get cluster id.
            . kmeansGen snd k
            . zip [0..] -- To keep track of index.
            . fmap S.toDenseListSV

-- | Get the cosine similarity between two rows using B2.
getSimilarityFromB2 :: B2 -> Int -> Int -> Double
getSimilarityFromB2 (B2 b2) i j =
    S.dot (S.extractRow b2 i) (S.extractRow b2 j)
        / (S.norm2 (S.extractRow b2 i) * S.norm2 (S.extractRow b2 j))

-- | Returns the eigenvector with the second smallest eigenvalue and on of the
-- symmetric normalized Laplacian L. Computes real symmetric part of L, so
-- ensure the input is real and symmetric. Diagonal should be 0s for adjacency
-- matrix. Uses I + Lnorm instead of I - Lnorm to find second largest singular
-- value instead of second smallest for Lnorm.
spectralNorm :: Int -> AdjacencyMatrix -> [S.SpVector Double]
spectralNorm e mat = secondLeft e lNorm
  where
    lNorm    = i S.^+^ (S.transpose invRootD S.#~# (mat S.#~# invRootD))
    invRootD = S.diagonalSM
             . S.vr
             . fmap ((\x -> if x == 0 then x else x ** (- 1 / 2)) . sum)
             . S.toRowsL
             $ mat
    i        = S.eye . S.nrows $ mat

-- | Returns the eigenvector with the second smallest eigenvalue and on of the
-- symmetric normalized Laplacian L. Computes real symmetric part of L, so
-- ensure the input is real and symmetric. Diagonal should be 0s for adjacency
-- matrix. Clusters the eigenvector using kmeans into k groups.
spectralClusterKNorm :: Int -> Int -> AdjacencyMatrix -> LabelVector
spectralClusterKNorm e k = kmeansVec k . spectralNorm e

-- | Returns the eigenvector with the second smallest eigenvalue of the
-- symmetric normalized Laplacian L. Computes real symmetric part of L, so
-- ensure the input is real and symmetric. Diagonal should be 0s for adjacency
-- matrix. Clusters the eigenvector by sign.
spectralClusterNorm :: AdjacencyMatrix -> LabelVector
spectralClusterNorm =
    fmap (bool 0 1 . (>= 0)) . foldl1' S.concatSV . spectralNorm 1
