static unsigned mul_step(unsigned A, unsigned B, unsigned M) {
    unsigned R = 0;
    unsigned n = M;
    while (n) {
        if (B&1) R += A;
        if (R&1) R += M;
        R >>= 1;
        B >>= 1;
        n >>= 1;
    }
    if (R >= M) R -= M;
    return R;
}

unsigned mul(unsigned A, unsigned B, unsigned M, unsigned B1) {
    return mul_step(mul_step(A, B, M), B1, M);
}