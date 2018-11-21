module MultiPol
  ( Polynomial(..)
  , Monomial(..)
  , toListOfMonomials
  , simplifiedListOfMonomials
  , fromListOfMonomials
  , toCanonicalForm
  , (^*^)
  , derivPoly
  , monom
  , evalPoly
  , polytest )
  where
import           Data.List
import           Data.Function
import qualified Data.Sequence as S
import           Data.Sequence (Seq, elemIndexL, (!?), adjust', findIndexL)
import           Data.Foldable (toList)

data Polynomial = Zero
                | M Monomial
                | Polynomial :+: Polynomial
                | Polynomial :*: Polynomial
                deriving(Show, Eq)

data Monomial = Monomial {
        coefficient :: Double,
        powers      :: Seq Int
    }
    deriving(Show, Eq)

zero :: Int -> Polynomial
zero n = M Monomial { coefficient = 0, powers = S.replicate n 0 }

derivMonomial :: Monomial -> Seq Int -> Polynomial
derivMonomial mono vars = let pows = powers mono in
  case sum vars of
    0 -> M mono
    1 -> let (Just i) = elemIndexL 1 vars in
         if pows !? i == Just 0 || coefficient mono == 0
           then Zero -- zero (S.length pows)
           else let (Just p) = pows !? i in
                M Monomial {coefficient = coefficient mono * fromIntegral p
                          , powers = adjust' (subtract 1) i pows
          }
    _ -> let (Just i) = findIndexL (>0) vars in
         let vars' = adjust' (subtract 1) i vars in
         let (Just p) = pows !? i in
         if pows !? i == Just 0 || coefficient mono == 0
           then Zero -- zero (S.length pows)
           else derivMonomial
                Monomial { coefficient = coefficient mono * fromIntegral p
                         , powers = adjust' (subtract 1) i pows }
                vars'

derivMonomial' :: Monomial -> [Int] -> Polynomial
derivMonomial' mono vars = derivMonomial mono (S.fromList vars)

multMonomial :: Monomial -> Monomial -> Monomial
multMonomial (Monomial ca powsa) (Monomial cb powsb) =
  Monomial (ca*cb) (S.zipWith (+) powsa powsb)

-- | polynomial to list of monomials
toListOfMonomials :: Polynomial -> [Monomial]
toListOfMonomials p = case p of
  M monomial -> if coefficient monomial == 0 then [] else [monomial]
  a :+: b -> toListOfMonomials a ++ toListOfMonomials b
  a :*: b -> [multMonomial monoa monob | monoa <- toListOfMonomials a,
                                         monob <- toListOfMonomials b]

-- | polynomial to list of monomials, grouping the monomials with same powers
simplifiedListOfMonomials :: Polynomial -> [Monomial]
simplifiedListOfMonomials p = map (foldl1 addMonomials) groups
  where
    groups = groupBy ((==) `on` powers)
             (sortBy (compare `on` powers) (toListOfMonomials p))
    -- cantorPairing :: Int -> Int -> Int
    -- cantorPairing k1 k2 = (k1+k2)*(k1+k2+1) + 2*k2
    -- cantorPairing' :: Seq Int -> Int
    -- cantorPairing' = foldl1 cantorPairing
    addMonomials :: Monomial -> Monomial -> Monomial
    addMonomials monoa monob = Monomial {
                             coefficient = coefficient monoa + coefficient monob
                           , powers = powers monoa
                           }

(.:+:.) :: Polynomial -> Polynomial -> Polynomial
(.:+:.) a b
  | a==Zero = b
  | b==Zero = a
  | otherwise = a :+: b

(.:*:.) :: Polynomial -> Polynomial -> Polynomial
(.:*:.) a b = if a == Zero || b == Zero
                 then Zero
                 else a :*: b

-- | differentiation
derivPoly :: Polynomial -> [Int] -> Polynomial
derivPoly pol vars = case pol of
  Zero    -> Zero
  M mono  -> derivMonomial' mono vars
  a :+: b -> derivPoly a vars .:+:. derivPoly b vars
  a :*: b -> (derivPoly a vars .:*:. b) .:+:. (a .:*:. derivPoly b vars)

-- | build a polynomial from a list of monomials
fromListOfMonomials :: [Monomial] -> Polynomial
fromListOfMonomials ms = if null ms
                            then Zero
                            else foldl1 (:+:) (map M ms)

-- | canonical form of a polynomial (sum of monomials)
toCanonicalForm :: Polynomial -> Polynomial
toCanonicalForm = fromListOfMonomials . simplifiedListOfMonomials

-- | scale polynomial by a scalar
(^*^) :: Double -> Polynomial -> Polynomial
(^*^) lambda pol = if lambda == 0
  then Zero
  else case pol of
    Zero -> Zero
    M monomial -> M (scaleMonomial monomial)
    a :+: b -> if a /= Zero && b /= Zero
      then (^*^) lambda a :+: (^*^) lambda b
      else if a == Zero
        then (^*^) lambda b
        else (^*^) lambda a
    a :*: b -> if a == Zero || b == Zero
      then Zero
      else (^*^) lambda a :*: b
  where
    scaleMonomial monomial = Monomial {
                                    coefficient = lambda * coefficient monomial
                                  , powers = powers monomial
                             }
-- | convenient built of a monomial
monom :: Double -> [Int] -> Monomial
monom coef pows = Monomial coef (S.fromList pows)

-- xpol :: Polynomial
-- xpol = M Monomial {coefficient = 1, powers = (1,0,0)}
--
-- ypol :: Polynomial
-- ypol = M Monomial {coefficient = 1, powers = (0,1,0)}

evalMonomial :: [Double] -> Monomial -> Double
evalMonomial xyz monomial =
  coefficient monomial * product (zipWith (^) xyz pows)
  where
    pows = toList (powers monomial)

evalPol :: Polynomial -> [Double] -> Double
evalPol pol xyz = sum (map (evalMonomial xyz) monomials)
  where
    monomials = toListOfMonomials pol

-- | evaluates a polynomial
evalPoly :: Polynomial -> [Double]-> Double
evalPoly pol xyz = case pol of
  Zero -> 0
  M mono -> evalMonomial xyz mono
  a :+: b -> evalPoly a xyz + evalPoly b xyz
  a :*: b -> evalPoly a xyz * evalPoly b xyz

-- evalFromListOfMonomials :: [Monomial] -> (Double, Double, Double) -> Double
-- evalFromListOfMonomials monomials xyz = sum (map (evalMonomial xyz) monomials)

polytest :: Polynomial
polytest = (M (monom 2 [3,1,1]) :+: M (monom 1 [2,0,0])) :*: M (monom 5 [1,1,1])
