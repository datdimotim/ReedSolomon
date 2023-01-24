package org.dimotim.reedsolomon;

import java.util.Arrays;
import java.util.OptionalInt;
import java.util.Random;
import java.util.function.IntBinaryOperator;
import java.util.function.IntPredicate;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class Main {
    public static void main(String[] args) {
        if(true){
            while (true) {
                long st = System.currentTimeMillis();
                for (int i = 0; i < 10; i++) {
                    check(97, 63);
                }
                System.out.println("--time-- " + (System.currentTimeMillis() - st));
            }
        } else {
            int[] dec = new ReedSolomon(new Field(97), 63).decode(new int[]{75, 28, 79, 59, 57, 61, 92, 66, 16, 1, 17, 40, 10, 76, 29, 14, 96, 84, 46, 0, 29, 27, 73, 83, 22, 48, 46, 29, 22, 53, 7, 53, 34, 68, 17, 7, 11, 34, 10, 66, 43, 4, 48, 9, 71, 56, 8, 34, 46, 79, 48, 20, 73, 85, 40, 46, 50, 65, 33, 76, 87, 20, 51, 18, 81, 87, 57, 94, 61, 29, 92, 8, 17, 83, 24, 71, 46, 47, 39, 62, 58, 81, 1, 24, 37, 17, 42, 31, 76, 83, 69, 47, 52, 43, 93, 76, 84});
            System.out.println(Arrays.equals(dec, new int[]{75, 60, 30, 69, 7, 25, 1, 24, 31, 22, 89, 36, 26, 90, 0, 84, 56, 73, 83, 32, 27, 23, 19, 63, 92, 29, 43, 53, 42, 81, 72, 18, 79, 68, 6, 67, 65, 59, 4, 19, 46, 89, 3, 74, 3, 27, 12, 71, 60, 4, 91, 19, 59, 11, 67, 60, 42, 52, 53, 91, 73, 31, 58}));
        }
    }

    private static void check(final int q, final int k) {
        final ReedSolomon reedSolomon = new ReedSolomon(new Field(q), k);
        final int e = (q - k)/2;

        //System.out.println("field size, code word length, q=" + q);
        //System.out.println("message length, k=" + k);
        //System.out.println("max errors, e=" + e);

        Random random = new Random();

        int[] src = IntStream.range(0, k).map(i -> random.nextInt(q)).toArray();
        //System.out.println("src word: " + Arrays.toString(src));

        int[] encoded = reedSolomon.encode(src);


        int[] sent = encoded.clone();
        IntStream.range(0, random.nextInt(e + 1))
                .map(i -> random.nextInt(q))
                .forEach(w -> sent[w] = random.nextInt(q));

        //System.out.println("encoded: " + Arrays.toString(encoded));
        //System.out.println("sent:    " + Arrays.toString(sent));

        int[] decoded = reedSolomon.decode(sent);
        //System.out.println("decoded: "+Arrays.toString(decoded));

        boolean success = Arrays.equals(src, decoded);
        System.out.println("success: "+ success);

        if(!success) {
            System.err.println("src: "+ Arrays.toString(src));
            System.err.println("sent: "+ Arrays.toString(sent));
            throw new RuntimeException("not success");
        }
    }
}

/**
 * @param f q = n - поле Fq, размер кодового слова; должно быть простым
 * @param k - количество информационных символов в кодовом слове
 * k = 2 * e, где e - максимальное количество ошибок в кодовом слове
 * k < n
 */
record ReedSolomon(Field f, int k, Poly poly) {
    ReedSolomon(Field f, int k) {
        this(f, k, new Poly(f));
        if (k >= f.q()) throw new IllegalArgumentException("k > q, k="+k + " q="+ f.q());
    }

    int q() {
        return f.q();
    }

    int[] encode(int[] src) {
        if (src.length != k) {
            throw new IllegalArgumentException("src length != k");
        }
        return IntStream.range(0, q())
                .map(x -> poly.at(x, src))
                .toArray();
    }

    int[] decode(int[] src) {
        if (src.length != q()) throw new IllegalArgumentException("src length != q");
        int[] pp = poly.interpolate(IntStream.range(0, q()).toArray(), src);

        int[] b = IntStream.range(0, q())
                .map(z -> f.mul(f.pow(z, (q() - k())/2), poly.at(z, pp)))
                .toArray();

        int[][] m = IntStream.range(0, q())
                .mapToObj(z -> IntStream.concat(
                        IntStream.rangeClosed(0, (q() + k())/2 - 1).map(j -> f.pow(z, j)),
                        IntStream.rangeClosed(0, (q() - k())/2 -1).map(j -> f.minus(0, f.mul(f.pow(z, j),poly.at(z, pp))))
                ).toArray())
                .toArray(int[][]::new);

        int[] xs = solve(m, b);

        int[] q = Arrays.stream(xs).limit((q()+k)/2).toArray();
        int[] r = IntStream.concat(Arrays.stream(xs).skip(q.length), IntStream.of(1)).toArray();

        int[] res = poly.divide(q, r);
        return IntStream.concat(Arrays.stream(res), IntStream.range(0, k - res.length).map(i -> 0)).toArray();
    }

    int[] solve(int[][] m, int[] b) {
        if (m[0].length == 1) {
            if (m[0][0] == 0) {
                if (b[0] != 0) throw new IllegalArgumentException("unsolvable system");
                return new int[]{0};
            } else {
                return new int[]{f.div(b[0], m[0][0])};
            }
        }

        if (m[0][0] == 0) {
            OptionalInt idx = IntStream.range(0, m.length).filter(i -> m[i][0] != 0).findFirst();
            if (idx.isEmpty()) {
                int[][] mSub = Arrays.stream(m).map(mi -> Arrays.stream(mi).skip(1).toArray())
                        .toArray(int[][]::new);
                return IntStream.concat(IntStream.of(0), Arrays.stream(solve(mSub, b))).toArray();
            } else {
                int[][] mm = swapM(0, idx.getAsInt(), m);
                int[] bb = swapB(0, idx.getAsInt(), b);
                return solve(mm, bb);
            }
        } else {
            int[][] mm = IntStream.range(1, m.length)
                    .mapToObj(i -> IntStream.range(1, m[i].length)
                            .map(j -> f.minus(m[i][j], f.mul(m[0][j], f.div(m[i][0], m[0][0]))))
                            .toArray())
                    .toArray(int[][]::new);

            int[] bb = IntStream.range(1, m.length)
                    .map(i -> f.minus(b[i], f.mul(b[0], f.div(m[i][0], m[0][0]))))
                    .toArray();

            int[] ss = solve(mm, bb);
            int x0 = f.div(IntStream.concat(
                        IntStream.of(b[0]),
                        IntStream.range(0, ss.length).map(i -> f.minus(0, f.mul(m[0][1 + i], ss[i])))
                    )
                    .reduce(f::plus).getAsInt(), m[0][0]);

            return IntStream.concat(IntStream.of(x0), Arrays.stream(ss)).toArray();
        }
    }

    void printMatr(int[][] m){
        System.out.println("matr:");
        Arrays.stream(m)
                .forEach(mi -> System.out.println(Arrays.stream(mi).mapToObj(mij -> mij + "\t").collect(Collectors.joining())));
        System.out.println("\n\n");
    }


    int[][] swapM(int i, int j, int[][] m) {
        m = m.clone();
        int[] mi = m[i];
        int[] mj = m[j];
        m[i] = mj;
        m[j] = mi;
        return m;
    }

    int[] swapB(int i, int j, int[] b) {
        b = b.clone();
        int bi = b[i];
        int bj = b[j];
        b[i] = bj;
        b[j] = bi;
        return b;
    }
}

/**
 *
 * @param f поле
 */
record Poly(Field f) {
    /**
     *
     * @param x точка
     * @param as коэффициенты начиная с младшей степени
     * @return значение в точке x
     */
    int at(int x, int[] as) {
        int v = 0;
        for (int i = 0; i < as.length; i++) {
            v = f.plus(v, f.mul(as[i],f.pow(x, i)));
        }
        return v;
    }

    int[] interpolate(int[] xs, int[] ys) {
        /*
                      (x  - x2)(x  - x3)...(x  - xn)
             f(x) =   ------------------------------ * y1 + ...
                      (x1 - x2)(x1 - x3)...(x1 - xn)

         */
        int n = xs.length;
        return IntStream.range(0, n)
                .mapToObj(i -> {
                    int xi = xs[i];
                    int zn = IntStream.range(0, n)
                            .filter(j -> j != i)
                            .map(j -> {
                                int xj = xs[j];
                                return  f.minus(xi, xj);
                            })
                            .reduce(f::mul).getAsInt();

                    int[] ch = IntStream.range(0, n)
                            .filter(j -> j != i)
                            .mapToObj(j -> {
                                int xj = xs[j];
                                return new int[]{f.minus(0, xj), 1};
                            }).reduce(this::mul).get();


                    int c = f.div(ys[i], zn);

                    return mulC(ch, c);
                }).reduce(this::plus).get();
    }

    int[] divide(int[] as, int[] bs) {
        if (as.length == 0) return new int[0];

       int gr = f.div(as[as.length-1], bs[bs.length-1]);

       int shift = as.length - bs.length;

       int[] asNext = Arrays.stream(minus(as, mulC(mulXn(bs, shift), gr))).limit(as.length - 1).toArray();

       return plus(divide(asNext, bs), mulXn(new int[]{gr}, shift));
    }

    int[] mul(int[] as, int[] bs) {
        if (bs.length == 1) {
            return mulC(as, bs[0]);
        }
        return plus(mulC(as, bs[0]), mulX(mul(as, Arrays.stream(bs).skip(1).toArray())));
    }

    int[] plusMinus(int[] as, int[] bs, IntBinaryOperator op){
        int[] rs = IntStream.range(0, Math.max(as.length, bs.length))
                .map(i -> op.applyAsInt((as.length > i ? as[i] : 0), (bs.length > i ? bs[i] : 0)))
                .toArray();

        int l = IntStream.range(0, rs.length)
                .map(i -> rs.length - 1 - i)
                .filter(i -> rs[i] != 0)
                .findFirst()
                .orElse(-1);

        return Arrays.stream(rs).limit(l + 1).toArray();
    }

    int[] minus(int[] as, int[] bs){
        return plusMinus(as, bs, f::minus);
    }

    int[] plus(int[] as, int[] bs){
        return plusMinus(as, bs, f::plus);
    }

    int[] mulX(int[] as) {
        return mulXn(as, 1);
    }

    int[] mulXn(int[] as, int n) {
        return IntStream.concat(IntStream.range(0, n).map(i -> 0), Arrays.stream(as)).toArray();
    }

    int[] mulC(int[] as, int c) {
        return Arrays.stream(as).map(a -> f.mul(a, c)).toArray();
    }
}

/**
 * поле Fq, q должно быть простым
 */
record Field (int q, int[] invs) {
    Field(int q) {
        this(q, getInvs(q));
        if(IntStream.range(2, q - 1).anyMatch(d -> q % d == 0)) {
            throw new IllegalArgumentException("q is not prime: " + q);
        }
    }

    private static int[] getInvs(int q) {
        return IntStream.concat(
                IntStream.of(0),
                IntStream.range(1, q)
                        .map(n -> IntStream.range(1, q)
                                .filter(i -> 1 == (i * n) % q)
                                .findFirst().getAsInt())
        ).toArray();
    }

    void assertVal(int v) {
        assertVal(v, x -> true);
    }

    void assertValAndNotNull(int v) {
        assertVal(v, x -> x != 0);
    }

    void assertVal(int v, IntPredicate p) {
        if (v < 0 || v >= q || !p.test(v)) {
            throw new IllegalArgumentException("invalid value: "+ v);
        }
    }

    int plus(int a, int b) {
        assertVal(a);
        assertVal(b);
        return (a + b) % q;
    }

    int minus(int a, int b) {
        assertVal(a);
        assertVal(b);
        return (q + a - b) % q;
    }

    int mul(int a, int b) {
        assertVal(a);
        assertVal(b);
        return (a * b) % q;
    }

    int div(int a, int b) {
        assertVal(a);
        assertValAndNotNull(b);
        return mul(a, invs[b]);
    }

    int pow(int a, int n) {
        assertVal(a);
        if (n < 0) {
            throw new IllegalArgumentException("n < 0, n = " + n);
        }

        if (a == 0) {
            return n != 0 ? 0 : 1;
        }

        return powH(a, n % (q - 1));  // a^phi(q) `mod` q = 1
    }

    private int powH(int a, int n) {
        if (n == 0) return 1;
        int r = powH(a, n / 2);
        int r2 =  mul(r, r);
        if (n % 2 == 0) {
            return r2;
        } else {
            return mul(a, r2);
        }
    }
}