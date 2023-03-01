package ckks;

import java.math.BigInteger;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Random;

public class Context {
    final static boolean debug = false;

    int slots;
    int topLevel;
    double defaultScale;

    // Roots of the cyclotomic polynomial
    Complex cyclotomicRoots[];

    // All generated primes, used to check that we didn't took some prime twice
    long[] allPrimes;
    int lastPrimeIdx = 0;

    // RNS primes
    long[] primes;

    // Level specific data
    LevelData[] levelsData;

    // Roots of unity modulo each prime
    long[][] primesRootsOfUnity;

    // VandermondeLU of each prime, used for encoding
    long[][][] primesVandermondeLUs;

    // in [i][j] the inverse of the j'th prime modulo the i'th prime
    long[][] otherPrimesInv;

    // temp primes. Used in relinearization
    long[] tempPrimes;

    // in [i] product of temp primes modulo the i'th regular prime
    long[] tempPrimesProdsMod;

    // in [i] inverse of product of temp primes modulo the i'th regular prime
    long[] tempPrimesProdsInvMod;

    // in [i][j] product of temp primes besides the i'th temp prime, modulo the j'th
    // regular prime
    long[][] otherTempPrimesProdsMods;

    // in [i] inverse of product of temp primes besides the i'th temp prime, modulo
    // the i'th temp prime
    long[] otherTempPrimesProdsInvsMod;

    /**
     * Initializes an empty context to be loaded with deserialize
     */
    public Context() {
    }

    /**
     * Initializes a context containing basic elements required for encoding,
     * encrypting, etc.
     * 
     * @param slots               Number of slots in each ciphertext. Must be a
     *                            power of 2.
     * @param multiplications     Number of multiplications the can be applied to a
     *                            ciphertext or a plaintext.
     * @param integerPrecision    Maximal number of bits (including the sign bit) of
     *                            the integer part of every number in a ciphertext
     *                            or a plaintext.
     * @param fractionalPrecision
     * @throws Exception
     */
    public Context(int slots, int multiplications, int integerPrecision, int fractionalPrecision) {
        if (!isPowerOfTwo(slots))
            throw new IllegalArgumentException("Number of slots must be power of 2");

        if (integerPrecision < 1)
            throw new IllegalArgumentException("integerPrecision must be at least 1");

        if (integerPrecision + fractionalPrecision > 62)
            throw new IllegalArgumentException("integerPrecision + fractionalPrecision can be at most 62");

        this.slots = slots;
        this.topLevel = multiplications;
        this.defaultScale = Math.pow(2, fractionalPrecision);

        if (debug)
            System.out.println("initCyclotomicRoots()\n");
        initCyclotomicRoots();

        if (debug)
            System.out.println("initPrimes()\n");
        initPrimes(multiplications, integerPrecision, fractionalPrecision);
        initPrimesData();
        initTempPrimes(fractionalPrecision);
        initTempPrimesData();

        if (debug)
            System.out.println("initPrimesRootsOfUnity()\n");
        initPrimesRootsOfUnity();

        if (debug)
            System.out.println("initPrimesVandermondeLUs()\n");
        initPrimesVandermondeLUs();
    }

    public void validateLevelDataExists(int level) {
        if (levelsData[level] == null)
            levelsData[level] = new LevelData(primes, tempPrimes, level);
    }

    private boolean isPowerOfTwo(int x) {
        // see:
        // https://stackoverflow.com/questions/600293/how-to-check-if-a-number-is-a-power-of-2
        return (x & (x - 1)) == 0;
    }

    private void initCyclotomicRoots() {
        cyclotomicRoots = new Complex[slots];

        for (int i = 0; i < slots; i++) {
            double theta = (Math.PI * (2 * i + 1)) / (2 * slots);
            cyclotomicRoots[i] = new Complex(Math.cos(theta), Math.sin(theta));
        }
    }

    private void initPrimes(int multiplications, int integerPrecision, int fractionalPrecision) {
        allPrimes = new long[2 * (multiplications + 1)];
        primes = new long[multiplications + 1];

        Random rnd = new Random();
        // We search for primes that are 1 bit larger than the wanted precision, so that
        // the range of valid values for encryption (determined by the precision
        // parameters) will be contained inside the modulos range. Also, notice that
        // JAVA doesn't have unsigned types, and so we limit our width to 63 bits max.
        primes[0] = genPrime(integerPrecision + fractionalPrecision + 1, rnd, slots * 2);
        for (int i = 1; i < primes.length; i++)
            primes[i] = genPrime(fractionalPrecision + 1, rnd, slots * 2);
    }

    private void initPrimesData() {
        if (debug) {
            System.out.println("Primes= " + Arrays.toString(primes) + '\n');

            BigInteger modulus = BigInteger.ONE;
            for (int i = 0; i < primes.length; i++)
                modulus = modulus.multiply(BigInteger.valueOf(primes[i]));

            System.out.println("Modulus= " + modulus + '\n');
        }

        levelsData = new LevelData[primes.length];

        otherPrimesInv = new long[primes.length][primes.length];
        for (int primeIdx = 0; primeIdx < primes.length; primeIdx++)
            for (int otherPrimeIdx = 0; otherPrimeIdx < primes.length; otherPrimeIdx++)
                if (primeIdx != otherPrimeIdx)
                    otherPrimesInv[primeIdx][otherPrimeIdx] = Maths.modInv(primes[otherPrimeIdx], primes[primeIdx]);
    }

    private void initTempPrimes(int fractionalPrecision) {
        tempPrimes = new long[primes.length];

        Random rnd = new Random();
        for (int i = 0; i < tempPrimes.length; i++)
            tempPrimes[i] = genPrime(fractionalPrecision, rnd, slots * 2);
    }

    private void initTempPrimesData() {

        if (debug) {
            System.out.println("Temp Primes= " + Arrays.toString(tempPrimes) + '\n');

            BigInteger modulus = BigInteger.ONE;
            for (int i = 0; i < tempPrimes.length; i++)
                modulus = modulus.multiply(BigInteger.valueOf(tempPrimes[i]));

            System.out.println("Temp Modulus= " + modulus + '\n');
        }

        tempPrimesProdsMod = new long[primes.length];
        for (int primeIdx = 0; primeIdx < primes.length; primeIdx++) {
            tempPrimesProdsMod[primeIdx] = 1;
            for (int tmpPrimeIdx = 0; tmpPrimeIdx < tempPrimes.length; tmpPrimeIdx++)
                tempPrimesProdsMod[primeIdx] = Maths.modMult(tempPrimesProdsMod[primeIdx], tempPrimes[tmpPrimeIdx],
                        primes[primeIdx]);
        }

        if (debug) {
            System.out.println("tempPrimesProdsMod= " + Arrays.toString(tempPrimesProdsMod) + '\n');
        }

        tempPrimesProdsInvMod = new long[primes.length];
        for (int primeIdx = 0; primeIdx < primes.length; primeIdx++)
            tempPrimesProdsInvMod[primeIdx] = Maths.modInv(tempPrimesProdsMod[primeIdx], primes[primeIdx]);

        if (debug) {
            System.out.println("tempPrimesProdsInvMod= " + Arrays.toString(tempPrimesProdsInvMod) + '\n');
        }

        otherTempPrimesProdsMods = new long[tempPrimes.length][primes.length];
        for (int primeIdx = 0; primeIdx < primes.length; primeIdx++) {
            for (int tempPrimeIdx = 0; tempPrimeIdx < tempPrimes.length; tempPrimeIdx++) {
                otherTempPrimesProdsMods[tempPrimeIdx][primeIdx] = 1;
                for (int i = 0; i < tempPrimes.length; i++) {
                    if (i == tempPrimeIdx)
                        continue;
                    otherTempPrimesProdsMods[tempPrimeIdx][primeIdx] = Maths
                            .modMult(otherTempPrimesProdsMods[tempPrimeIdx][primeIdx], tempPrimes[i], primes[primeIdx]);
                }
            }
        }

        if (debug) {
            System.out.println("otherTempPrimesProdsMods= " + Arrays.deepToString(otherTempPrimesProdsMods) + '\n');
        }

        otherTempPrimesProdsInvsMod = new long[tempPrimes.length];
        for (int tempPrimeIdx = 0; tempPrimeIdx < tempPrimes.length; tempPrimeIdx++) {
            otherTempPrimesProdsInvsMod[tempPrimeIdx] = 1;
            for (int i = 0; i < tempPrimes.length; i++) {
                if (i == tempPrimeIdx)
                    continue;
                otherTempPrimesProdsInvsMod[tempPrimeIdx] = Maths.modMult(otherTempPrimesProdsInvsMod[tempPrimeIdx],
                        tempPrimes[i], tempPrimes[tempPrimeIdx]);
            }
            otherTempPrimesProdsInvsMod[tempPrimeIdx] = Maths.modInv(otherTempPrimesProdsInvsMod[tempPrimeIdx],
                    tempPrimes[tempPrimeIdx]);
        }

        if (debug) {
            System.out.println("otherTempPrimesProdsInvsMod= " + Arrays.toString(otherTempPrimesProdsInvsMod) + '\n');
        }

    }

    private long genPrime(int bits, Random rnd, int N) {
        // for each prime p, (p-1) needs to be divisible by 2N (slots * 4)
        // so that the cyclotomic polynomial can be factorized completely modulo p.
        // see:
        // https://en.wikipedia.org/wiki/Cyclotomic_polynomial#Cyclotomic_polynomials_over_a_finite_field_and_over_the_p-adic_integers
        BigInteger M = BigInteger.valueOf(2 * N);

        long res = 0;
        boolean found = false;
        int tries = 0;

        while (!found && tries < 100) {
            tries++;

            BigInteger num = BigInteger.probablePrime(bits, rnd);

            if (!num.subtract(BigInteger.ONE).mod(M).equals(BigInteger.ZERO))
                continue;

            res = num.longValue();

            found = true;

            for (int i = 0; i < lastPrimeIdx; i++)
                if (allPrimes[i] == res)
                    found = false;
        }

        if (!found)
            throw new IllegalStateException(
                    "Unable to find enough primes for given parameters. Try to increase fractional precision.");

        allPrimes[lastPrimeIdx++] = res;

        return res;
    }

    // see:
    // https://math.stackexchange.com/questions/158344/how-to-find-the-solutions-for-the-n-th-root-of-unity-in-modular-arithmetic
    private void initPrimesRootsOfUnity() {
        int N = slots * 2;

        primesRootsOfUnity = new long[primes.length][N];

        Random rand = new Random();
        long generator = 0, res;
        for (int primeIdx = 0; primeIdx < primes.length; primeIdx++) {
            HashSet<Long> generators = new HashSet<Long>();
            res = 1;

            // if res == 1 then we didn't got a generator
            while (res == 1) {
                generator = Maths.mod(rand.nextLong(), primes[primeIdx] - 2) + 1;
                if (generators.contains(generator))
                    continue;

                generators.add(generator);
                res = Maths.modPow(generator, (primes[primeIdx] - 1) / 2, primes[primeIdx]);
            }

            for (int rootIdx = 0; rootIdx < N; rootIdx++)
                primesRootsOfUnity[primeIdx][rootIdx] = Maths.modPow(generator,
                        ((primes[primeIdx] - 1) / (2 * N)) * (2 * rootIdx + 1), primes[primeIdx]);

            if (debug) {
                System.out.println("(x^" + N + " + 1) roots modulo " + primes[primeIdx] + ": "
                        + Arrays.toString(primesRootsOfUnity[primeIdx]) + '\n');
            }

        }
    }

    private void initPrimesVandermondeLUs() {
        int N = slots * 2;
        primesVandermondeLUs = new long[primes.length][N][N];

        for (int primeIdx = 0; primeIdx < primes.length; primeIdx++) {
            // init vandermondes
            for (int rootIdx = 0; rootIdx < N; rootIdx++) {
                long val = 1;
                for (int exp = 0; exp < N; exp++) {
                    primesVandermondeLUs[primeIdx][rootIdx][exp] = val;
                    val = Maths.modMult(val, primesRootsOfUnity[primeIdx][rootIdx], primes[primeIdx]);
                }
            }

            if (debug) {
                System.out.println("vandermonde matrix for prime " + primes[primeIdx]);
                System.out.println(Arrays.deepToString(primesVandermondeLUs[primeIdx]));
                System.out.println();
            }

            // LU decomposition
            // based on code from https://math.nist.gov/javanumerics/jama/
            long[] LUrowi;
            long[] LUcolj = new long[N];

            // Outer loop.
            for (int j = 0; j < N; j++) {

                // Make a copy of the j-th column to localize references.
                for (int i = 0; i < N; i++)
                    LUcolj[i] = primesVandermondeLUs[primeIdx][i][j];

                // Apply previous transformations.
                for (int i = 0; i < N; i++) {
                    LUrowi = primesVandermondeLUs[primeIdx][i];

                    int kmax = Math.min(i, j);
                    long a[] = new long[kmax];
                    long b[] = new long[kmax];

                    for (int k = 0; k < kmax; k++) {
                        a[k] = LUrowi[k];
                        b[k] = LUcolj[k];
                    }

                    long s = 0;
                    for (int k = 0; k < kmax; k++) {
                        s += Maths.modMult(LUrowi[k], LUcolj[k], primes[primeIdx]);
                        s = Maths.mod(s, primes[primeIdx]);
                    }

                    LUcolj[i] -= s;
                    LUcolj[i] = Maths.mod(LUcolj[i], primes[primeIdx]);
                    LUrowi[j] = LUcolj[i];
                }

                // Compute multipliers.
                if (j < N & primesVandermondeLUs[primeIdx][j][j] != 0) {
                    for (int i = j + 1; i < N; i++)
                        primesVandermondeLUs[primeIdx][i][j] = Maths.modDiv(primesVandermondeLUs[primeIdx][i][j],
                                primesVandermondeLUs[primeIdx][j][j], primes[primeIdx]);

                }
            }

            if (debug) {
                System.out.println("LU decomposition of vandermonde matrix for prime " + primes[primeIdx]);
                System.out.println(Arrays.deepToString(primesVandermondeLUs[primeIdx]));
                System.out.println();
            }
        }
    }

    public long[] serialize() {
        long[] res = new long[5 + primes.length + tempPrimes.length
                + primesRootsOfUnity.length * primesRootsOfUnity[0].length];

        int idx = 0;

        res[idx++] = slots;
        res[idx++] = topLevel;
        res[idx++] = (long) defaultScale;
        res[idx++] = primes.length;
        res[idx++] = tempPrimes.length;

        for (int i = 0; i < primes.length; i++)
            res[idx++] = primes[i];

        for (int i = 0; i < tempPrimes.length; i++)
            res[idx++] = tempPrimes[i];

        for (int i = 0; i < primes.length; i++)
            for (int j = 0; j < slots * 2; j++)
                res[idx++] = primesRootsOfUnity[i][j];

        return res;
    }

    public void deserialize(long[] serialization) {
        int idx = 0;

        slots = (int) serialization[idx++];
        topLevel = (int) serialization[idx++];
        defaultScale = (double) serialization[idx++];
        primes = new long[(int) serialization[idx++]];
        tempPrimes = new long[(int) serialization[idx++]];

        initCyclotomicRoots();

        for (int i = 0; i < primes.length; i++)
            primes[i] = serialization[idx++];

        initPrimesData();

        for (int i = 0; i < tempPrimes.length; i++)
            tempPrimes[i] = serialization[idx++];

        initTempPrimesData();

        primesRootsOfUnity = new long[primes.length][slots * 2];

        for (int i = 0; i < primes.length; i++)
            for (int j = 0; j < slots * 2; j++)
                primesRootsOfUnity[i][j] = serialization[idx++];

        initPrimesVandermondeLUs();
    }

    public int getNumSlots() {
        return slots;
    }
}
