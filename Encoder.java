package ckks;

import java.lang.Math;
import java.util.Arrays;

public class Encoder {
    final boolean debug = false;

    private Context context;

    // this is also the basis for the ring of vector of polynomial evaluations
    // (sigma(R))
    private Complex[][] vandermondeTransposed;

    private Complex[][] vandermondeLU;
    private int[] piv;

    public Encoder(Context context) {
        this.context = context;

        calcVandermondeTransposed(context.slots);

        if (debug) {
            System.out.println("Vandermonde transposed:");
            System.out.println(Arrays.deepToString(vandermondeTransposed));
            System.out.println();
        }

        calcVandermondeLU();

        if (debug) {
            System.out.println("VandermondeLU:");
            System.out.println(Arrays.deepToString(vandermondeLU));
            System.out.println();
        }
    }

    public void encode(Complex[] src, Plaintext res) {
        if (src.length > context.slots)
            throw new IllegalArgumentException("Vec size should be at most " + context.slots);

        Complex[] vec = expandVector(src);
        if (debug) {
            System.out.println("Expended vector:");
            System.out.println(Arrays.toString(vec));
            System.out.println();
        }

        scaleVector(vec);
        if (debug) {
            System.out.println("Scaled vector:");
            System.out.println(Arrays.toString(vec));
            System.out.println();
        }

        discretizeVector(vec);
        if (debug) {
            System.out.println("Discretized vector:");
            System.out.println(Arrays.toString(vec));
            System.out.println();
            System.out.println("Sanity check:");
            Complex[] sanityCheck = Complex.multVectorByScalar(vec, 1.0 / context.defaultScale);
            System.out.println(Arrays.toString(sanityCheck));
            System.out.println();
        }

        long[] coeffs = canonicalEmbeddingInverse(vec);
        if (debug) {
            System.out.println("Polynomial coefficients:");
            System.out.println(Arrays.toString(coeffs));
            System.out.println();
            System.out.println("Sanity check:");
            double[] sanityCheckCoeffs = new double[coeffs.length];
            for (int i = 0; i < sanityCheckCoeffs.length; i++)
                sanityCheckCoeffs[i] = (double) coeffs[i] / context.defaultScale;
            System.out.println("scaled coeffs:");
            System.out.println(Arrays.toString(sanityCheckCoeffs));
            Complex[] sanityCheck = new Complex[context.slots];
            for (int i = 0; i < context.slots; i++) {
                double theta = (Math.PI * (2 * i + 1)) / (2 * context.slots);
                Complex x = new Complex(Math.cos(theta), Math.sin(theta));
                Complex x_copy = new Complex(x);
                sanityCheck[i] = new Complex(sanityCheckCoeffs[0]);
                for (int j = 1; j < sanityCheckCoeffs.length; j++) {
                    sanityCheck[i].add_inplace(x_copy.mult(new Complex(sanityCheckCoeffs[j])));
                    x_copy.mult_inplace(x);
                }
            }
            System.out.println(Arrays.toString(sanityCheck));
            System.out.println();
        }

        long[][] rnsCoeffs = Maths.rns(context, coeffs);
        Maths.ntt_inplace(context, rnsCoeffs, context.primes.length - 1);

        Polynomial poly = new Polynomial(context, rnsCoeffs);

        if (debug) {
            System.out.println("polynomial after RNS and NTT");
            poly.debugPrint();
            System.out.println();
        }

        res.init(poly, context.topLevel);
    }

    public Complex[] decode(Plaintext p) {
        if (debug) {
            System.out.println("Plaintext level= " + p.getLevel());
            System.out.println("Plaintext scale= " + p.getScale());
            System.out.println();
        }

        long[][] crt = p.getM().getCrt();

        if (debug) {
            System.out.println("Plaintext polynomial double CRT reps.");
            System.out.println(Arrays.deepToString(crt));
            System.out.println();
        }

        long[][] rnsCoeffs = new long[context.primes.length][context.slots * 2];
        for (int i = 0; i < context.primes.length; i++)
            for (int j = 0; j < context.slots * 2; j++)
                rnsCoeffs[i][j] = crt[i][j];

        Maths.nttInverse_inplace(context, rnsCoeffs, p.getLevel());

        if (debug) {
            System.out.println("Plaintext polynomial after NTT inverse");
            System.out.println(Arrays.deepToString(rnsCoeffs));
            System.out.println();
        }

        double[] coeffs = Maths.rnsInverse(context, rnsCoeffs, p.getLevel());

        if (debug) {
            System.out.println("Plaintext polynomial after RNS inverse");
            System.out.println(Arrays.toString(coeffs));
            System.out.println();
        }

        for (int i = 0; i < coeffs.length; i++)
            coeffs[i] /= p.getScale();

        if (debug) {
            System.out.println("Plaintext polynomial after scaling");
            System.out.println(Arrays.toString(coeffs));
            System.out.println();
        }

        // TODO: use FFT
        // evaluate on roots of the cyclotomic polynomial
        Complex[] res = new Complex[context.slots];

        for (int i = 0; i < context.slots; i++) {
            Complex x = context.cyclotomicRoots[i];

            Complex x_copy = new Complex(x);

            res[i] = new Complex(coeffs[0]);

            for (int j = 1; j < coeffs.length; j++) {
                res[i].add_inplace(x_copy.mult(new Complex(coeffs[j])));
                x_copy.mult_inplace(x);
            }
        }

        return res;
    }

    private void calcVandermondeTransposed(int slots) {
        vandermondeTransposed = new Complex[slots * 2][];

        for (int exp = 0; exp < slots * 2; exp++) {
            Complex[] row = new Complex[slots * 2];

            for (int i = 0; i < slots * 2; i++) {
                double theta = (Math.PI * (2 * i + 1) * exp) / (2 * slots);
                row[i] = new Complex(Math.cos(theta), Math.sin(theta));
            }

            vandermondeTransposed[exp] = row;
        }
    }

    // based on code from https://math.nist.gov/javanumerics/jama/
    private void calcVandermondeLU() {
        int N = context.slots * 2;

        vandermondeLU = new Complex[N][N];

        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
                vandermondeLU[j][i] = new Complex(vandermondeTransposed[i][j]);

        piv = new int[N];

        for (int i = 0; i < N; i++) {
            piv[i] = i;
        }

        int pivsign = 1;

        Complex[] LUrowi;
        Complex[] LUcolj = new Complex[N];

        // Outer loop.
        for (int j = 0; j < N; j++) {

            // Make a copy of the j-th column to localize references.
            for (int i = 0; i < N; i++)
                LUcolj[i] = new Complex(vandermondeLU[i][j]);

            // Apply previous transformations.
            for (int i = 0; i < N; i++) {
                LUrowi = vandermondeLU[i];

                int kmax = Math.min(i, j);
                Complex a[] = new Complex[kmax];
                Complex b[] = new Complex[kmax];

                for (int k = 0; k < kmax; k++) {
                    a[k] = LUrowi[k];
                    b[k] = LUcolj[k];
                }

                Complex s = new Complex(0, 0);
                for (int k = 0; k < kmax; k++) {
                    s.add_inplace(LUrowi[k].mult(LUcolj[k]));
                }
                LUcolj[i].sub_inplace(s);
                LUrowi[j] = new Complex(LUcolj[i]);
            }

            // Find pivot and exchange if necessary.
            int p = j;
            for (int i = j + 1; i < N; i++) {
                if (LUcolj[i].norm() > LUcolj[p].norm()) {
                    p = i;
                }
            }
            if (p != j) {
                for (int k = 0; k < N; k++) {
                    Complex t = new Complex(vandermondeLU[p][k]);
                    vandermondeLU[p][k] = new Complex(vandermondeLU[j][k]);
                    vandermondeLU[j][k] = new Complex(t);
                }
                int k = piv[p];
                piv[p] = piv[j];
                piv[j] = k;
                pivsign = -pivsign;
            }

            // Compute multipliers.
            if (j < N & !vandermondeLU[j][j].equal(new Complex(0, 0))) {
                for (int i = j + 1; i < N; i++) {
                    vandermondeLU[i][j].divide_inplace(vandermondeLU[j][j]);
                }
            }
        }
    }

    private Complex[] expandVector(Complex[] vec) {
        Complex[] expendedVec = new Complex[context.slots * 2];

        for (int i = 0; i < context.slots; i++)
            if (i < vec.length)
                expendedVec[i] = new Complex(vec[i]);
            else
                expendedVec[i] = new Complex(0, 0);

        for (int i = context.slots - 1; i >= 0; i--)
            expendedVec[context.slots * 2 - 1 - i] = expendedVec[i].conj();

        return expendedVec;
    }

    private void scaleVector(Complex[] vec) {
        for (int i = 0; i < vec.length; i++)
            vec[i].mult_inplace(new Complex(context.defaultScale));
    }

    private void discretizeVector(Complex[] vec) {
        double[] coordinates = calcCoordinates(vec);
        if (debug) {
            System.out.println("Coordinates of vector in sigma(R)'s basis:");
            System.out.println(Arrays.toString(coordinates));

            System.out.println("Sanity check:");
            Complex[] sanityCheck = Complex.multVectorByScalar(vandermondeTransposed[0], coordinates[0]);
            for (int i = 1; i < vandermondeTransposed.length; i++)
                Complex.addVectors_inplace(sanityCheck,
                        Complex.multVectorByScalar(vandermondeTransposed[i], coordinates[i]));
            System.out.println(Arrays.toString(sanityCheck));
            System.out.println();
        }

        long[] roundedCoordinates = roundCoordinatesRandomly(coordinates);
        if (debug) {
            System.out.println("Coordinates rounded randomly:");
            System.out.println(Arrays.toString(roundedCoordinates));
            System.out.println();
        }

        // Generate discretized vector from basis and rounded coordinates
        Complex.multVectorByScalar_inplace(vec, 0);
        for (int i = 0; i < vandermondeTransposed.length; i++)
            Complex.addVectors_inplace(vec,
                    Complex.multVectorByScalar(vandermondeTransposed[i], roundedCoordinates[i]));
    }

    private double[] calcCoordinates(Complex[] vec) {
        double[] coordinates = new double[vec.length];

        for (int i = 0; i < vec.length; i++) {
            double hermProd = Complex.hermitianProduct(vec, vandermondeTransposed[i]).re;

            double normSquared = Complex.hermitianProduct(vandermondeTransposed[i], vandermondeTransposed[i]).re;

            coordinates[i] = hermProd / normSquared;
        }

        return coordinates;
    }

    private long[] roundCoordinatesRandomly(double[] coordinates) {
        long[] res = new long[coordinates.length];

        for (int i = 0; i < coordinates.length; i++) {
            double elem = coordinates[i];
            long elemFloored = (long) Math.floor(elem);

            double distanceFromFloor = elem - elemFloored;

            // uniform distribution over [0,1]
            double rand = Math.random();

            // see
            // https://stackoverflow.com/questions/40183948/how-to-generate-random-number-based-on-probability-in-java
            if (distanceFromFloor < rand)
                res[i] = elemFloored;
            else
                res[i] = elemFloored + 1;
        }

        return res;
    }

    private long[] canonicalEmbeddingInverse(Complex[] vec) {
        long[] coeffs = new long[context.slots * 2];

        // replace with matrix multiplication
        Complex[] solution = solve(vec);

        for (int i = 0; i < coeffs.length; i++)
            coeffs[i] = Math.round(solution[i].re);

        return coeffs;
    }

    // based on code from https://math.nist.gov/javanumerics/jama/
    private Complex[] solve(Complex[] vec) {
        int N = context.slots * 2;

        // Copy right hand side with pivoting
        Complex[] X = new Complex[vec.length];

        for (int i = 0; i < piv.length; i++)
            X[i] = new Complex(vec[piv[i]]);

        // Solve L*Y = B(piv,:)
        for (int k = 0; k < N; k++) {
            for (int i = k + 1; i < N; i++) {
                X[i].sub_inplace(X[k].mult(vandermondeLU[i][k]));
            }
        }
        // Solve U*X = Y;
        for (int k = N - 1; k >= 0; k--) {
            X[k].divide_inplace(vandermondeLU[k][k]);

            for (int i = 0; i < k; i++)
                X[i].sub_inplace(X[k].mult(vandermondeLU[i][k]));
        }

        return X;
    }
}
