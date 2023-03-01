package ckks;

import java.math.BigInteger;
import java.util.Arrays;

public class Maths {
   final static boolean debug = false;

   public static long[][] ntt(Context context, long[][] coeffs, int level) {
      long[][] res = new long[coeffs.length][coeffs[0].length];
      for (int i = 0; i < coeffs.length; i++)
         for (int j = 0; j < coeffs[0].length; j++)
            res[i][j] = coeffs[i][j];

      ntt_inplace(context, res, level);
      return res;
   }

   public static void ntt_inplace(Context context, long[][] coeffs, int level) {
      // TODO: change this to do ntt only on used levels

      if (debug) {
         System.out.println("Maths.ntt_inplace:");
         System.out.println("coeffs:");
         System.out.println(Arrays.deepToString(coeffs));
      }

      int N = context.slots * 2;

      long variable, res;
      long[] primeCoeffs = new long[N];

      // TODO: use FFT
      for (int primeIdx = 0; primeIdx <= level; primeIdx++) {
         for (int coeffIdx = 0; coeffIdx < N; coeffIdx++)
            primeCoeffs[coeffIdx] = coeffs[primeIdx][coeffIdx];

         for (int rootIdx = 0; rootIdx < N; rootIdx++) {
            variable = context.primesRootsOfUnity[primeIdx][rootIdx];
            res = primeCoeffs[0];

            for (int exp = 1; exp < primeCoeffs.length; exp++) {
               res += modMult(variable, primeCoeffs[exp], context.primes[primeIdx]);
               res = mod(res, context.primes[primeIdx]);

               if (exp < primeCoeffs.length - 1)
                  variable = modMult(variable, context.primesRootsOfUnity[primeIdx][rootIdx], context.primes[primeIdx]);
            }

            coeffs[primeIdx][rootIdx] = res;
         }
      }

      if (debug) {
         System.out.println("ntt:");
         System.out.println(Arrays.deepToString(coeffs));
         System.out.println();
      }
   }

   public static long[][] nttInverse(Context context, long[][] ntt, int level) {
      long[][] res = new long[ntt.length][ntt[0].length];
      for (int i = 0; i < ntt.length; i++)
         for (int j = 0; j < ntt[0].length; j++)
            res[i][j] = ntt[i][j];

      nttInverse_inplace(context, res, level);
      return res;
   }

   public static void nttInverse_inplace(Context context, long[][] ntt, int level) {
      if (debug) {
         System.out.println("Maths.nttInverse_inplace");
         System.out.println("ntt:");
         System.out.println(Arrays.deepToString(ntt));
      }

      int N = context.slots * 2;

      // based on code from https://math.nist.gov/javanumerics/jama/
      for (int primeIdx = 0; primeIdx <= level; primeIdx++) {
         // Solve L*Y = B
         for (int k = 0; k < N; k++) {
            for (int i = k + 1; i < N; i++) {
               ntt[primeIdx][i] -= modMult(ntt[primeIdx][k], context.primesVandermondeLUs[primeIdx][i][k],
                     context.primes[primeIdx]);
               ntt[primeIdx][i] = mod(ntt[primeIdx][i], context.primes[primeIdx]);
            }
         }

         if (debug) {
            System.out.println("Coeffs modulo " + context.primes[primeIdx] + " after solving L*Y=B:");
            System.out.println(Arrays.toString(ntt[primeIdx]));
         }

         // Solve U*X = Y;
         for (int k = N - 1; k >= 0; k--) {
            ntt[primeIdx][k] = modDiv(ntt[primeIdx][k], context.primesVandermondeLUs[primeIdx][k][k],
                  context.primes[primeIdx]);

            for (int i = 0; i < k; i++) {
               ntt[primeIdx][i] -= modMult(ntt[primeIdx][k], context.primesVandermondeLUs[primeIdx][i][k],
                     context.primes[primeIdx]);
               ntt[primeIdx][i] = mod(ntt[primeIdx][i], context.primes[primeIdx]);
            }
         }

         if (debug) {
            System.out.println("Coeffs modulo " + context.primes[primeIdx] + " after solving U*X = Y:");
            System.out.println(Arrays.toString(ntt[primeIdx]));
         }
      }

      if (debug) {
         System.out.println();
      }
   }

   public static long[][] rns(Context context, long[] coeffs) {
      if (debug) {
         System.out.println("Maths.rns");
         System.out.println("coeffs:");
         System.out.println(Arrays.toString(coeffs));
      }

      int N = context.slots * 2;

      long[][] res = new long[context.primes.length][N];

      for (int primeIdx = 0; primeIdx < context.primes.length; primeIdx++)
         for (int coeffIdx = 0; coeffIdx < N; coeffIdx++)
            res[primeIdx][coeffIdx] = mod(coeffs[coeffIdx], context.primes[primeIdx]);

      if (debug) {
         System.out.println("rns coeffs:");
         System.out.println(Arrays.deepToString(res));
         System.out.println();
      }

      return res;
   }

   public static double[] rnsInverse(Context context, long[][] coeffs, int level) {
      if (debug) {
         System.out.println("Maths.rnsInverse");
         System.out.println("rns coeefs:");
         System.out.println(Arrays.deepToString(coeffs));
      }

      context.validateLevelDataExists(level);
      LevelData levelData = context.levelsData[level];

      int N = context.slots * 2;

      double[] res = new double[N];

      for (int i = 0; i < N; i++) {
         BigInteger coeff = BigInteger.valueOf(0);

         // TODO: need to calculate those for every level
         for (int primeIdx = 0; primeIdx <= level; primeIdx++)
            coeff = coeff.add(levelData.otherPrimesProds[primeIdx]
                  .multiply(
                        levelData.otherPrimesProdsInvs[primeIdx].multiply(BigInteger.valueOf(coeffs[primeIdx][i]))));

         coeff = coeff.mod(levelData.primesProd);
         if (coeff.compareTo(levelData.halfPrimesProd) == 1) {
            if (debug)
               System.out.println(
                     "coeff #" + i + " (" + coeff.toString() + ") is bigger than "
                           + levelData.halfPrimesProd.toString());

            coeff = coeff.subtract(levelData.primesProd);
         }

         res[i] = coeff.doubleValue();
      }

      if (debug) {
         System.out.println("final coeffs:");
         System.out.println(Arrays.toString(res));
         System.out.println();
      }

      return res;
   }

   public static BigInteger[] experimental_rnsInverse(Context context, long[][] coeffs, int level) {
      if (debug) {
         System.out.println("Maths.experimental_rnsInverse");
         System.out.println("rns coeefs:");
         System.out.println(Arrays.deepToString(coeffs));
      }

      context.validateLevelDataExists(level);
      LevelData levelData = context.levelsData[level];

      int N = context.slots * 2;

      BigInteger[] res = new BigInteger[N];

      for (int i = 0; i < N; i++) {
         BigInteger coeff = BigInteger.valueOf(0);

         // TODO: need to calculate those for every level
         for (int primeIdx = 0; primeIdx <= level; primeIdx++)
            coeff = coeff.add(levelData.otherPrimesProds[primeIdx]
                  .multiply(
                        levelData.otherPrimesProdsInvs[primeIdx].multiply(BigInteger.valueOf(coeffs[primeIdx][i]))));

         coeff = coeff.mod(levelData.primesProd);

         res[i] = coeff;
      }

      if (debug) {
         System.out.println("final coeffs:");
         System.out.println(Arrays.toString(res));
         System.out.println();
      }

      return res;
   }

   public static long[][] modUp(Context context, long[][] d2Coeffs, int level) {
      if (debug) {
         System.out.println("Maths.modUp");
         System.out.println("d2Coeffs:");
         System.out.println(Arrays.deepToString(d2Coeffs));
      }

      context.validateLevelDataExists(level);
      LevelData levelData = context.levelsData[level];

      long[][] newCoeffs = new long[context.primes.length + context.tempPrimes.length][context.slots * 2];

      for (int i = 0; i < context.primes.length; i++)
         for (int j = 0; j < context.slots * 2; j++)
            newCoeffs[i][j] = d2Coeffs[i][j];

      // fast basis conversion
      for (int tempPrimeIdx = 0; tempPrimeIdx < context.tempPrimes.length; tempPrimeIdx++) {
         for (int coeffIdx = 0; coeffIdx < context.slots * 2; coeffIdx++) {
            long value = 0, tmp = 0;
            for (int primeIdx = 0; primeIdx <= level; primeIdx++) {
               tmp = Maths.modMult(d2Coeffs[primeIdx][coeffIdx], levelData.otherPrimesProdsInvsMod[primeIdx],
                     context.primes[primeIdx]);
               tmp = Maths.modMult(tmp, levelData.otherPrimesProdsMods[primeIdx][tempPrimeIdx],
                     context.tempPrimes[tempPrimeIdx]);
               value = Maths.mod(value + tmp, context.tempPrimes[tempPrimeIdx]);
            }

            newCoeffs[context.primes.length + tempPrimeIdx][coeffIdx] = value;
         }
      }

      if (debug) {
         System.out.println("newCoeffs:");
         System.out.println(Arrays.deepToString(newCoeffs));
         System.out.println();
      }

      return newCoeffs;
   }

   public static long[][] modDown(Context context, long[][] d2Coeffs, int level) {
      if (debug) {
         System.out.println("Maths.modDown");
         System.out.println("d2Coeffs:");
         System.out.println(Arrays.deepToString(d2Coeffs));
      }

      long[][] newCoeffs = new long[context.primes.length][context.slots * 2];

      // fast basis conversion
      for (int primeIdx = 0; primeIdx <= level; primeIdx++) {
         for (int coeffIdx = 0; coeffIdx < context.slots * 2; coeffIdx++) {
            long value = 0, tmp = 0;
            for (int tempPrimeIdx = 0; tempPrimeIdx < context.tempPrimes.length; tempPrimeIdx++) {
               tmp = Maths.modMult(d2Coeffs[context.primes.length + tempPrimeIdx][coeffIdx],
                     context.otherTempPrimesProdsInvsMod[tempPrimeIdx], context.tempPrimes[tempPrimeIdx]);
               tmp = Maths.modMult(tmp, context.otherTempPrimesProdsMods[tempPrimeIdx][primeIdx],
                     context.primes[primeIdx]);
               value = Maths.mod(value + tmp, context.primes[primeIdx]);
            }

            value = Maths.mod(d2Coeffs[primeIdx][coeffIdx] - value, context.primes[primeIdx]);

            value = Maths.modMult(value, context.tempPrimesProdsInvMod[primeIdx], context.primes[primeIdx]);

            newCoeffs[primeIdx][coeffIdx] = value;
         }
      }

      if (debug) {
         System.out.println("newCoeffs:");
         System.out.println(Arrays.deepToString(newCoeffs));
         System.out.println();
      }

      return newCoeffs;
   }

   // the polynomials are in coefficients form. we need to calculate their product
   // mod m but also mod (x^N+1)
   public static long[] multPolynominalsMod(long[] a, long[] b, long m) {
      if (debug) {
         System.out.println("Maths.multPolynominalsMod");
         System.out.println("a= " + Arrays.toString(a));
         System.out.println("b= " + Arrays.toString(b));
         System.out.println("m= " + m);
      }

      int N = a.length;

      long[] res = new long[N];

      long tmp = 0;
      int idx = 0;
      for (int i = 0; i < a.length; i++) {
         for (int j = 0; j < b.length; j++) {
            idx = (int) Maths.mod(i + j, N);
            tmp = Maths.modMult(a[i], b[j], m);
            if (i + j >= N)
               tmp = -tmp;
            res[idx] = Maths.mod(res[idx] + tmp, m);
         }
      }

      if (debug) {
         System.out.println(Arrays.toString(res));
         System.out.println();
      }

      return res;
   }

   public static long mod(long a, long m) {
      if (debug) {
         System.out.println("Maths.mod");
         System.out.println("Maths.mod with a= " + a + ", m= " + m);
      }

      long res = a % m;

      if (debug) {
         System.out.println("Maths.mod res= " + (res >= 0 ? res : res + m) + "\n");
      }

      return res >= 0 ? res : res + m;
   }

   public static long modInv(long a, long m) {
      if (debug) {
         System.out.println("Maths.modInv with a= " + a + ", m= " + m);
      }

      a = mod(a, m);

      long[] gcdRes = xgcd(a, m);

      if (gcdRes[0] != 1)
         throw new IllegalStateException("Number " + a + " is not invertible modulo " + m);

      long res = mod(gcdRes[1], m);

      if (debug) {
         System.out.println("Maths.modInv res= " + res + "\n");
      }

      return res;
   }

   // mod multiplication without overflow. see:
   // https://stackoverflow.com/questions/20971888/modular-multiplication-of-large-numbers-in-c
   public static long modMult(long a, long z, long m) {
      if (debug) {
         System.out.println("Maths.modMult with a= " + a + ", z=" + z + ", m= " + m);
      }

      a = mod(a, m);
      z = mod(z, m);

      long res;

      // if multiplying a and z cause overflow multiplyExact will throw
      try {
         res = mod(Math.multiplyExact(a, z), m);
      } catch (ArithmeticException e) {
         long q = m / a;
         long r = mod(m, a);
         res = a * mod(z, q) - modMult(r, (z / q), m);
         if (res < 0)
            res += m;
      }

      if (debug) {
         System.out.println("Maths.modMult res= " + res + "\n");
      }

      return res;

   }

   public static long modDiv(long a, long z, long m) {
      if (debug) {
         System.out.println("Maths.modDiv with a= " + a + ", z=" + z + ", m= " + m);
      }

      long inv = modInv(z, m);

      long res = modMult(a, inv, m);

      if (debug) {
         System.out.println("Maths.modDiv res= " + res + "\n");
      }

      return res;
   }

   public static long modPow(long b, long e, long m) {
      if (debug) {
         System.out.println("Maths.modPow with b= " + b + ", e=" + e + ", m= " + m);
      }

      b = mod(b, m);

      long x = 1;
      long curr_base = b;
      while (e > 0) {
         if (e % 2 == 1)
            x = modMult(x, curr_base, m);

         curr_base = modMult(curr_base, curr_base, m);
         e /= 2;
      }

      long res = mod(x, m);

      if (debug) {
         System.out.println("Maths.modPow res= " + res + "\n");
      }

      return res;
   }

   public static long[] ternaryDist(int N, double rho) {
      long[] res = new long[N];

      double rand;
      for (int i = 0; i < N; i++) {
         // uniform distribution over [0,1]
         rand = Math.random();

         // see
         // https://stackoverflow.com/questions/40183948/how-to-generate-random-number-based-on-probability-in-java
         if (rand < rho / 2)
            res[i] = -1;
         else if (rand < rho)
            res[i] = 1;
         else
            res[i] = 0;
      }

      return res;
   }

   private static long[] xgcd(long a, long b) {
      if (a == 0)
         return new long[] { b, 0, 1 };

      if (a == 1)
         return new long[] { 1, b + 1, -1 };

      long[] tmp = xgcd(mod(b, a), a);
      return new long[] { tmp[0], tmp[2] - tmp[1] * (b / a), tmp[1] };
   }
}
