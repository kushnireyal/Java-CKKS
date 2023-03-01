package ckks;

import java.math.BigInteger;
import java.util.Arrays;

public class LevelData {
    final static boolean debug = false;

    // Product of all prime in this level
    BigInteger primesProd;

    // Half the above value
    BigInteger halfPrimesProd;

    // in [i] the product of all primes except for the i'th prime
    BigInteger[] otherPrimesProds;

    // in [i] the inverse of the product of all primes except for the i'th prime,
    // modulo the i'th prime
    BigInteger[] otherPrimesProdsInvs;

    // in [i][j] product of regular primes besides the i'th regular prime, modulo
    // the j'th temp prime
    long[][] otherPrimesProdsMods;

    // in [i] inverse of product of regular primes besides the i'th regular prime,
    // modulo the i'th regular prime
    long[] otherPrimesProdsInvsMod;

    public LevelData(long[] primes, long[] tempPrimes, int level) {
        primesProd = BigInteger.valueOf(1);
        for (int primeIdx = 0; primeIdx <= level; primeIdx++)
            primesProd = primesProd.multiply(BigInteger.valueOf(primes[primeIdx]));
        halfPrimesProd = primesProd.divide(BigInteger.valueOf(2));

        otherPrimesProds = new BigInteger[primes.length];
        otherPrimesProdsInvs = new BigInteger[primes.length];
        for (int primeIdx = 0; primeIdx <= level; primeIdx++) {
            otherPrimesProds[primeIdx] = primesProd.divide(BigInteger.valueOf(primes[primeIdx]));
            otherPrimesProdsInvs[primeIdx] = otherPrimesProds[primeIdx].mod(BigInteger.valueOf(primes[primeIdx]))
                    .modInverse(BigInteger.valueOf(primes[primeIdx]));
        }

        otherPrimesProdsMods = new long[primes.length][tempPrimes.length];
        for (int tempPrimeIdx = 0; tempPrimeIdx < tempPrimes.length; tempPrimeIdx++) {
            BigInteger m = BigInteger.valueOf(tempPrimes[tempPrimeIdx]);
            for (int primeIdx = 0; primeIdx <= level; primeIdx++) {
                otherPrimesProdsMods[primeIdx][tempPrimeIdx] = otherPrimesProds[primeIdx].mod(m).longValue();
            }
        }

        if (debug) {
            System.out.println("otherPrimesProdsMods=");
            System.out.println(Arrays.deepToString(otherPrimesProdsMods) + '\n');
        }

        otherPrimesProdsInvsMod = new long[primes.length];
        for (int primeIdx = 0; primeIdx <= level; primeIdx++)
            otherPrimesProdsInvsMod[primeIdx] = otherPrimesProdsInvs[primeIdx].mod(BigInteger.valueOf(primes[primeIdx]))
                    .longValue();

        if (debug) {
            System.out.println("otherPrimesProdsInvsMod=");
            System.out.println(Arrays.toString(otherPrimesProdsInvsMod) + '\n');
        }
    }
}
