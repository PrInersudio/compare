static unsigned bit_length(unsigned num) {
    unsigned count = 0;
    while (num) {
        count++;
        num >>= 1;
    }
    return count;
}

static unsigned mul_step(unsigned A, unsigned B, unsigned M, unsigned n) {
    unsigned R = 0;
    for (unsigned i = 0; i < n; ++i) {
        if (B&(1<<i)) R += A;
        if (R&1) R+=M;
        R >>= 1;
    }
    if (R >= M) R -= M;
    return R;
}

unsigned mul(unsigned A, unsigned B, unsigned M) {
    unsigned n = bit_length(M);
    return mul_step(mul_step(A, B, M, n), (1 << (2*n)) % M, M, n);
}