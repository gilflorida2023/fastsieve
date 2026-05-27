#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>

typedef unsigned __int128 u128;

#define FILENAME "primes_state.bin"
#define WHEEL_MOD 210
#define WHEEL_SIZE 48
#define WHEEL_BYTES 48

static const uint64_t wheel_residues[WHEEL_SIZE] = {
    1, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
    53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
    101, 103, 107, 109, 113, 121, 127, 131,
    137, 139, 143, 149, 151, 157, 163, 167, 169, 173,
    179, 181, 187, 191, 193, 197, 199, 209
};

static uint8_t residue_to_bit[WHEEL_MOD];

static void build_lut(void) {
    for (int i = 0; i < WHEEL_MOD; i++)
        residue_to_bit[i] = 0xFF;
    for (int i = 0; i < WHEEL_SIZE; i++)
        residue_to_bit[wheel_residues[i]] = i;
}

typedef struct {
    u128 prime;
    u128 next;
} PrimeEntry;

static PrimeEntry *base = NULL;
static uint64_t base_count = 0;
static uint64_t base_cap = 0;
static u128 base_sieved = 0;

static uint64_t total_count = 0;

static void print_u128(u128 n) {
    if (n == 0) { putchar('0'); return; }
    char buf[48];
    int i = 47;
    buf[i] = '\0';
    while (n > 0) {
        buf[--i] = '0' + (char)(n % 10);
        n /= 10;
    }
    fputs(buf + i, stdout);
}

static void print_u128_f(FILE *fp, u128 n) {
    if (n == 0) { fputc('0', fp); return; }
    char buf[48];
    int i = 47;
    buf[i] = '\0';
    while (n > 0) {
        buf[--i] = '0' + (char)(n % 10);
        n /= 10;
    }
    fputs(buf + i, fp);
}

static u128 sqrt_u128(u128 n) {
    if (n < 2) return n;
    u128 lo = 1, hi = (u128)1 << 64;
    if (hi > n) hi = n;
    while (lo < hi) {
        u128 mid = lo + (hi - lo) / 2;
        if (mid <= n / mid && (mid + 1) > n / (mid + 1))
            return mid;
        if (mid > n / mid)
            hi = mid;
        else
            lo = mid + 1;
    }
    return lo;
}

static int parse_u128(const char *s, u128 *out) {
    *out = 0;
    if (!s || !*s) return 0;
    while (*s) {
        if (*s < '0' || *s > '9') return 0;
        *out = *out * 10 + (u128)(*s - '0');
        s++;
    }
    return 1;
}

static int write_entry(FILE *fp, PrimeEntry *e) {
    uint64_t lo, hi;
    lo = (uint64_t)e->prime;
    hi = (uint64_t)(e->prime >> 64);
    if (fwrite(&lo, 8, 1, fp) != 1) return 0;
    if (fwrite(&hi, 8, 1, fp) != 1) return 0;
    lo = (uint64_t)e->next;
    hi = (uint64_t)(e->next >> 64);
    if (fwrite(&lo, 8, 1, fp) != 1) return 0;
    if (fwrite(&hi, 8, 1, fp) != 1) return 0;
    return 1;
}

static int read_entry(FILE *fp, PrimeEntry *e) {
    uint64_t lo, hi;
    if (fread(&lo, 8, 1, fp) != 1) return 0;
    if (fread(&hi, 8, 1, fp) != 1) return 0;
    e->prime = ((u128)hi << 64) | lo;
    if (fread(&lo, 8, 1, fp) != 1) return 0;
    if (fread(&hi, 8, 1, fp) != 1) return 0;
    e->next = ((u128)hi << 64) | lo;
    return 1;
}

static void add_base_prime(u128 p) {
    if (base_count >= base_cap) {
        base_cap = base_cap ? base_cap * 2 : 65536;
        base = realloc(base, base_cap * sizeof(PrimeEntry));
    }
    base[base_count].prime = p;
    base[base_count].next = 0;
    base_count++;
}

static void extend_base_primes(u128 limit) {
    if (limit <= base_sieved) return;
    if (limit > (u128)1 << 62) limit = (u128)1 << 62;

    if (base_count == 0) {
        add_base_prime(2);
        add_base_prime(3);
        add_base_prime(5);
        add_base_prime(7);
        base_sieved = 7;
        if (limit <= 7) return;
    }

    while (base_sieved < limit) {
        u128 next_target;
        if (base_sieved >= ((u128)1 << 62))
            next_target = limit;
        else
            next_target = base_sieved * base_sieved;
        if (next_target > limit || next_target <= base_sieved)
            next_target = limit;
        if (next_target <= base_sieved) break;

        u128 start = base_sieved + 1;
        u128 end = next_target;
        if (start > end) break;

        uint64_t first_block = (uint64_t)(start / WHEEL_MOD);
        uint64_t last_block = (uint64_t)(end / WHEEL_MOD);
        uint64_t num_blocks = last_block - first_block + 1;

        uint8_t *buf = malloc(num_blocks * WHEEL_BYTES);
        if (!buf) { fprintf(stderr, "Out of memory\n"); exit(1); }
        memset(buf, 0xFF, num_blocks * WHEEL_BYTES);

        for (uint64_t bi = 0; bi < base_count; bi++) {
            u128 p = base[bi].prime;
            if (p < 11) continue;
            u128 fm = start;
            if (p * p > fm) fm = p * p;
            if (fm > end) continue;
            if (fm % p != 0)
                fm += p - fm % p;
            u128 m = fm;
            while (m <= end) {
                uint64_t r = (uint64_t)(m % WHEEL_MOD);
                uint8_t ri = residue_to_bit[r];
                if (ri < WHEEL_SIZE) {
                    uint64_t block = (uint64_t)(m / WHEEL_MOD) - first_block;
                    buf[block * WHEEL_BYTES + ri] = 0;
                }
                m += p;
            }
        }

        for (uint64_t bi = 0; bi < num_blocks; bi++) {
            for (int ri = 0; ri < WHEEL_SIZE; ri++) {
                if (!buf[bi * WHEEL_BYTES + ri]) continue;
                u128 n = (first_block + bi) * (u128)WHEEL_MOD + wheel_residues[ri];
                if (n < start || n > end) continue;

                add_base_prime(n);

                u128 start_m = n * n;
                if (start_m < start) start_m = start;
                for (u128 m = start_m; m <= end; m += n) {
                    uint64_t r = (uint64_t)(m % WHEEL_MOD);
                    uint8_t ri2 = residue_to_bit[r];
                    if (ri2 < WHEEL_SIZE) {
                        uint64_t b = (uint64_t)(m / WHEEL_MOD) - first_block;
                        if (b < num_blocks)
                            buf[b * WHEEL_BYTES + ri2] = 0;
                    }
                }
            }
        }

        free(buf);
        base_sieved = end;
    }
}

static void process_segment(FILE *state_fp, u128 seg_start, u128 seg_end,
                            FILE *output_fp) {
    if (seg_start > seg_end) return;

    uint64_t first_block = (uint64_t)(seg_start / WHEEL_MOD);
    uint64_t last_block = (uint64_t)(seg_end / WHEEL_MOD);
    uint64_t num_blocks = last_block - first_block + 1;

    uint8_t *buf = malloc(num_blocks * WHEEL_BYTES);
    if (!buf) { fprintf(stderr, "Out of memory\n"); exit(1); }
    memset(buf, 0xFF, num_blocks * WHEEL_BYTES);

    for (uint64_t i = 0; i < base_count; i++) {
        u128 p = base[i].prime;
        u128 fm = seg_start;
        if (p * p > fm) fm = p * p;
        if (fm > seg_end) continue;
        if (fm % p != 0)
            fm += p - fm % p;

        u128 m = fm;
        while (m <= seg_end) {
            uint64_t r = (uint64_t)(m % WHEEL_MOD);
            uint8_t ri = residue_to_bit[r];
            if (ri < WHEEL_SIZE) {
                uint64_t block = (uint64_t)(m / WHEEL_MOD) - first_block;
                buf[block * WHEEL_BYTES + ri] = 0;
            }
            m += p;
        }
    }

    for (uint64_t bi = 0; bi < num_blocks; bi++) {
        for (int ri = 0; ri < WHEEL_SIZE; ri++) {
            if (!buf[bi * WHEEL_BYTES + ri]) continue;
            u128 n = (first_block + bi) * (u128)WHEEL_MOD + wheel_residues[ri];
            if (n < seg_start || n > seg_end) continue;

            total_count++;

            PrimeEntry e;
            e.prime = n;
            e.next = n * 2;
            while (e.next <= seg_end)
                e.next += n;

            if (state_fp) {
                fseek(state_fp, 0, SEEK_END);
                write_entry(state_fp, &e);
            }

            if (output_fp) {
                print_u128_f(output_fp, n);
                fputc('\n', output_fp);
            }

            u128 start_m = n * n;
            if (start_m < seg_start) start_m = seg_start;
            if (start_m > seg_end) continue;
            for (u128 m = start_m; m <= seg_end; m += n) {
                uint64_t r = (uint64_t)(m % WHEEL_MOD);
                uint8_t ri2 = residue_to_bit[r];
                if (ri2 < WHEEL_SIZE) {
                    uint64_t b = (uint64_t)(m / WHEEL_MOD) - first_block;
                    if (b < num_blocks)
                        buf[b * WHEEL_BYTES + ri2] = 0;
                }
            }
        }
    }

    free(buf);
}

static u128 auto_opt_buffer(u128 target) {
    u128 buf = ((u128)32768 / WHEEL_BYTES) * WHEEL_MOD;
    if (buf < WHEEL_MOD) buf = WHEEL_MOD;
    if (buf > target) buf = target;
    return buf;
}

static void format_duration(double sec, char *buf, size_t sz) {
    if (sec < 60.0) {
        snprintf(buf, sz, "%.4f sec", sec);
    } else if (sec < 3600.0) {
        int m = (int)(sec / 60);
        double s = sec - m * 60;
        snprintf(buf, sz, "%dm %.1fs", m, s);
    } else if (sec < 86400.0) {
        int h = (int)(sec / 3600);
        int m = (int)((sec - h * 3600) / 60);
        int s = (int)(sec - h * 3600 - m * 60);
        snprintf(buf, sz, "%dh %dm %ds", h, m, s);
    } else {
        int d = (int)(sec / 86400);
        int h = (int)((sec - d * 86400) / 3600);
        int m = (int)((sec - d * 86400 - h * 3600) / 60);
        snprintf(buf, sz, "%dd %dh %dm", d, h, m);
    }
}

static void print_help(void) {
    printf("Usage: primes [buffer_size] [target] [-c] [-r] [-o file]\n");
    printf("  buffer_size  segment window size in natural numbers (optional, auto-opt)\n");
    printf("  target       sieve up to this number (required for sieve; omit for -r)\n");
    printf("  -c           count only (no state file, faster)\n");
    printf("  -r           report from existing primes_state.bin (no sieve)\n");
    printf("  -o file      write primes to file (use \"-\" for stdout)\n");
}

int main(int argc, char **argv) {
    build_lut();

    u128 target = 0;
    u128 buffer_size = 0;
    char *output_filename = NULL;
    int count_only = 0;
    int report_mode = 0;

    int num_count = 0;
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-c") == 0) {
            count_only = 1;
        } else if (strcmp(argv[i], "-r") == 0) {
            report_mode = 1;
        } else if (strcmp(argv[i], "-o") == 0) {
            if (i + 1 < argc) output_filename = argv[++i];
            else { fprintf(stderr, "-o requires a filename\n"); return 1; }
        } else if (argv[i][0] == '-') {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            print_help();
            return 1;
        } else {
            u128 val;
            if (!parse_u128(argv[i], &val)) {
                fprintf(stderr, "Invalid number: %s\n", argv[i]);
                return 1;
            }
            if (num_count == 0) buffer_size = val;
            else if (num_count == 1) target = val;
            else { fprintf(stderr, "Unexpected argument: %s\n", argv[i]); return 1; }
            num_count++;
        }
    }
    if (num_count == 1) { target = buffer_size; buffer_size = 0; }

    if (count_only && report_mode) {
        fprintf(stderr, "Cannot use -c and -r together\n");
        return 1;
    }

    if (report_mode) {
        FILE *rfp = fopen(FILENAME, "rb");
        if (!rfp) { perror("Failed to open " FILENAME); return 1; }
        fseek(rfp, 0, SEEK_END);
        long fsize = ftell(rfp);
        rewind(rfp);

        PrimeEntry e;
        u128 first = 0, last = 0;
        uint64_t nentries = 0;
        while (read_entry(rfp, &e)) {
            if (nentries == 0) first = e.prime;
            last = e.prime;
            nentries++;
        }

        printf("State file: %s (%ld bytes)\n", FILENAME, fsize);
        printf("Entries:    %llu\n", (unsigned long long)nentries);
        printf("First:      ");
        print_u128(first);
        printf("\nLast:       ");
        print_u128(last);
        printf("\n");

        if (output_filename) {
            FILE *ofp;
            if (strcmp(output_filename, "-") == 0) {
                ofp = stdout;
            } else {
                ofp = fopen(output_filename, "w");
                if (!ofp) { perror(output_filename); return 1; }
            }
            rewind(rfp);
            while (read_entry(rfp, &e)) {
                print_u128_f(ofp, e.prime);
                fputc('\n', ofp);
            }
            if (ofp != stdout) fclose(ofp);
        }
        fclose(rfp);
        return 0;
    }

    if (target == 0) {
        print_help();
        return 1;
    }

    if (buffer_size == 0)
        buffer_size = auto_opt_buffer(target);

    printf("Target: ");
    print_u128(target);
    printf("\nBuffer: ");
    print_u128(buffer_size);
    printf("\n");

    FILE *state_fp = NULL;
    if (!count_only) {
        state_fp = fopen(FILENAME, "wb+");
        if (!state_fp) {
            perror("Failed to open state file");
            return 1;
        }
    }

    FILE *output_fp = NULL;
    if (output_filename) {
        output_fp = fopen(output_filename, "w");
        if (!output_fp) {
            perror("Failed to open output file");
            fclose(state_fp);
            return 1;
        }
    }

    clock_t start_time = clock();

    u128 current = 2;
    u128 sqrt_cur = 2;
    uint64_t segs = 0;

    static const u128 small_primes[] = {2, 3, 5, 7};
    for (int i = 0; i < 4; i++) {
        if (small_primes[i] > target) continue;
        total_count++;
        if (state_fp) {
            PrimeEntry e = {small_primes[i], small_primes[i] * 2};
            write_entry(state_fp, &e);
        }
        if (output_fp) {
            print_u128_f(output_fp, small_primes[i]);
            fputc('\n', output_fp);
        }
    }

    while (current <= target) {
        u128 seg_end = current + buffer_size - 1;
        if (seg_end > target) seg_end = target;

        sqrt_cur = sqrt_u128(seg_end);
        if (sqrt_cur > base_sieved)
            extend_base_primes(sqrt_cur);

        process_segment(state_fp, current, seg_end, output_fp);

        current = seg_end + 1;
        segs++;
    }

    clock_t end_time = clock();
    double elapsed = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    printf("Found ");
    print_u128(total_count);
    char dur[64];
    format_duration(elapsed, dur, sizeof(dur));
    printf(" primes in %llu segments (%s)\n",
           (unsigned long long)segs, dur);

    if (state_fp) fclose(state_fp);
    if (output_fp) fclose(output_fp);

    free(base);
    return 0;
}
