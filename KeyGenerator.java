package ckks;

import java.util.Arrays;
import java.util.concurrent.ThreadLocalRandom;

public class KeyGenerator {
    final static boolean debug = false;

    private Context context;

    private SecretKey secretKey;

    private long[] sCoeffs;

    private PublicKeys publicKeys;

    public KeyGenerator(Context context) {
        this.context = context;

        genSecretKey();
        genPublicKeys();
        genRelinearizationKey();
    }

    public SecretKey getSecretKey() {
        return secretKey;
    }

    public PublicKeys getPublicKeys() {
        return publicKeys;
    }

    public static Polynomial ternaryDist(Context context, double rho) {
        int N = context.slots * 2;
        long[] coeffs = Maths.ternaryDist(N, rho);

        if (debug) {
            System.out.println("ternaryDist");
            System.out.println(Arrays.toString(coeffs));
        }

        long[][] rnsCoeffs = Maths.rns(context, coeffs);
        Maths.ntt_inplace(context, rnsCoeffs, context.primes.length - 1);

        return new Polynomial(context, rnsCoeffs);
    }

    private void genSecretKey() {
        sCoeffs = Maths.ternaryDist(context.slots * 2, 0.5);

        long[][] rnsCoeffs = Maths.rns(context, sCoeffs);
        Maths.ntt_inplace(context, rnsCoeffs, context.primes.length - 1);

        Polynomial s = new Polynomial(context, rnsCoeffs);

        if (debug) {
            System.out.println("s coeffs:");
            System.out.println(Arrays.toString(sCoeffs));
            System.out.println();
            System.out.println("s:");
            s.debugPrint();
            System.out.println();
        }

        secretKey = new SecretKey(s);
    }

    private void genPublicKeys() {
        // uniform distribution
        long[][] aCrt = new long[context.primes.length][context.slots * 2];
        for (int primeIdx = 0; primeIdx < context.primes.length; primeIdx++)
            for (int i = 0; i < context.slots * 2; i++)
                aCrt[primeIdx][i] = ThreadLocalRandom.current().nextLong(0, context.primes[primeIdx]);

        if (debug) {
            long[][] aCrtCpy = new long[context.primes.length][context.slots * 2];
            for (int primeIdx = 0; primeIdx < context.primes.length; primeIdx++)
                for (int i = 0; i < context.slots * 2; i++)
                    aCrtCpy[primeIdx][i] = aCrt[primeIdx][i];

            Maths.nttInverse_inplace(context, aCrtCpy, context.primes.length - 1);

            System.out.println("pubA as polynomial:");
            System.out.println(
                    Arrays.toString(Maths.experimental_rnsInverse(context, aCrtCpy, context.primes.length - 1)));
        }

        Polynomial a = new Polynomial(context, aCrt);
        if (debug) {
            System.out.println("\npubE as polynomial:");
        }
        Polynomial e = ternaryDist(context, 0.5);
        Polynomial b = e.sub(a.mult(secretKey.getS()));
        if (debug) {
            System.out.println("pubE:");
            e.debugPrint();
            System.out.println("\npubB:");
            b.debugPrint();
            System.out.println("\npubA:");
            a.debugPrint();
            System.out.println();
        }

        publicKeys = new PublicKeys(b, a);
    }

    private void genRelinearizationKey() {
        int N = context.slots * 2;

        long[][] a = new long[context.primes.length + context.tempPrimes.length][N];

        // uniform distribution
        for (int i = 0; i < N; i++) {
            for (int primeIdx = 0; primeIdx < context.primes.length; primeIdx++)
                a[primeIdx][i] = ThreadLocalRandom.current().nextLong(0, context.primes[primeIdx]);

            for (int tempPrimeIdx = 0; tempPrimeIdx < context.tempPrimes.length; tempPrimeIdx++)
                a[context.primes.length + tempPrimeIdx][i] = ThreadLocalRandom.current().nextLong(0,
                        context.tempPrimes[tempPrimeIdx]);
        }

        Polynomial s = secretKey.getS();
        long[][] sSquaredRns = Maths.nttInverse(context, s.mult(s).getCrt(),
                context.primes.length - 1);

        if (debug) {
            System.out.println("s^2 RNS:");
            System.out.println(Arrays.deepToString(sSquaredRns) + '\n');
        }

        for (int primeIdx = 0; primeIdx < context.primes.length; primeIdx++)
            for (int coeffIdx = 0; coeffIdx < N; coeffIdx++)
                sSquaredRns[primeIdx][coeffIdx] = Maths.modMult(sSquaredRns[primeIdx][coeffIdx],
                        context.tempPrimesProdsMod[primeIdx], context.primes[primeIdx]);

        if (debug) {
            System.out.println("P*s^2 RNS (where P is prod of temp primes):");
            System.out.println(Arrays.deepToString(sSquaredRns) + '\n');
        }

        long[] e = Maths.ternaryDist(N, 0.5);

        if (debug) {
            System.out.println("RelinKeyE:");
            System.out.println(Arrays.toString(e) + '\n');
        }

        long[][] minusA = new long[context.primes.length + context.tempPrimes.length][N];
        for (int i = 0; i < a.length; i++)
            for (int j = 0; j < a[0].length; j++)
                minusA[i][j] = -a[i][j];

        if (debug) {
            System.out.println("-a:");
            System.out.println(Arrays.deepToString(minusA) + '\n');
        }

        long[][] b = new long[context.primes.length + context.tempPrimes.length][N];

        for (int primeIdx = 0; primeIdx < context.primes.length; primeIdx++) {
            b[primeIdx] = Maths.multPolynominalsMod(minusA[primeIdx], sCoeffs, context.primes[primeIdx]);

            if (debug) {
                System.out.println("-a^s mod " + context.primes[primeIdx] + ":");
                System.out.println(Arrays.toString(b[primeIdx]) + '\n');
            }

            for (int i = 0; i < e.length; i++) {
                b[primeIdx][i] = Maths.mod(b[primeIdx][i] + e[i], context.primes[primeIdx]);
                b[primeIdx][i] = Maths.mod(b[primeIdx][i] + sSquaredRns[primeIdx][i], context.primes[primeIdx]);
            }
        }

        for (int tempPrimeIdx = 0; tempPrimeIdx < context.tempPrimes.length; tempPrimeIdx++) {
            b[context.primes.length + tempPrimeIdx] = Maths.multPolynominalsMod(
                    minusA[context.primes.length + tempPrimeIdx], sCoeffs, context.tempPrimes[tempPrimeIdx]);
            for (int i = 0; i < e.length; i++)
                b[context.primes.length + tempPrimeIdx][i] = Maths
                        .mod(b[context.primes.length + tempPrimeIdx][i] + e[i], context.tempPrimes[tempPrimeIdx]);
        }

        if (debug) {
            System.out.println("RelinKeyB:");
            System.out.println(Arrays.deepToString(b) + '\n');
            System.out.println("RelinKeyA:");
            System.out.println(Arrays.deepToString(a) + '\n');
        }

        publicKeys.setRelinKeys(b, a);
    }
}
