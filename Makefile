CFLAGS = -O3 -pipe -Wall -Wextra
LDFLAGS =
TARGET = fastsieve

.PHONY: all clean test bench-wheels psvh-gen psvh-verify

all: $(TARGET)

$(TARGET): fastsieve.c
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^

clean:
	rm -f $(TARGET) primes_state.bin primes_state.ckpt

# Accuracy tests: all wheels at 10³, 10⁴, 10⁶
test: $(TARGET)
	@echo "=== Accuracy: π(1000) = 168 ==="
	@for w in 2 6 30 210 2310; do \
		got=$$(./$(TARGET) --wheel $$w 1000 -c 2>&1 | sed -n 's/Found \([0-9]*\) primes.*/\1/p'); \
		[ "$$got" = "168" ] || { echo "FAIL: wheel=$$w π(1000)=$$got"; exit 1; }; \
		echo "  PASS: wheel=$$w π(1000)=168"; \
	done
	@echo "=== Accuracy: π(10000) = 1229 ==="
	@for w in 2 6 30 210 2310; do \
		got=$$(./$(TARGET) --wheel $$w 10000 -c 2>&1 | sed -n 's/Found \([0-9]*\) primes.*/\1/p'); \
		[ "$$got" = "1229" ] || { echo "FAIL: wheel=$$w π(10⁴)=$$got"; exit 1; }; \
		echo "  PASS: wheel=$$w π(10⁴)=1229"; \
	done
	@echo "=== Accuracy: π(1000000) = 78498 ==="
	@for w in 2 6 30 210 2310; do \
		got=$$(./$(TARGET) --wheel $$w 1000000 -c 2>&1 | sed -n 's/Found \([0-9]*\) primes.*/\1/p'); \
		[ "$$got" = "78498" ] || { echo "FAIL: wheel=$$w π(10⁶)=$$got"; exit 1; }; \
		echo "  PASS: wheel=$$w π(10⁶)=78498"; \
	done
	@echo "=== Resume test (wheel 210) ==="
	@rm -f primes_state.bin primes_state.ckpt
	@./$(TARGET) --wheel 210 -s 1000000 -c >/dev/null 2>&1
	@got=$$(./$(TARGET) --wheel 210 -R 2000000 -c 2>&1 | sed -n 's/Found \([0-9]*\) primes.*/\1/p'); \
	[ "$$got" = "148933" ] || { echo "FAIL: resume π(2M)=$$got, expected 148933"; exit 1; }; \
	echo "  PASS: resume π(2M)=148933"
	@echo "=== Resume wheel mismatch test ==="
	@rm -f primes_state.bin primes_state.ckpt
	@./$(TARGET) --wheel 210 -s 1000000 -c >/dev/null 2>&1
	@output=$$(./$(TARGET) --wheel 30 -R 2000000 -c 2>&1 || true); \
	echo "$$output" | grep -q "wheel mismatch" || { echo "FAIL: wheel mismatch not detected"; exit 1; }; \
	echo "  PASS: wheel mismatch correctly errors"

# Benchmark: median of 3 runs at 10⁹, conditional 10¹² if <30s
bench-wheels: $(TARGET)
	@echo "=== Benchmark @ 10⁹ (median of 3 runs) ==="
	@for w in 30 210 2310; do \
		echo "  Wheel $$w:"; \
		times=""; \
		for i in 1 2 3; do \
			t=$$(/usr/bin/time -f "%e" ./$(TARGET) --wheel $$w 1000000000 -c 2>&1 >/dev/null); \
			times="$$times $$t"; \
		done; \
		med=$$(for t in $$times; do echo $$t; done | sort -n | head -2 | tail -1); \
		echo "    median: $$med s"; \
		echo "$$w $$med" >> /tmp/bench_times.$$; \
	done
	@echo ""
	@echo "=== Conditional Benchmark @ 10¹² (if 10⁹ < 30s) ==="
	@while read w med; do \
		if echo "$$med < 30" | bc -l | grep -q 1; then \
			echo "  Wheel $$w: 10⁹=$$med s < 30s → running 10¹²..."; \
			t=$$(/usr/bin/time -f "%e" timeout 1800 ./$(TARGET) --wheel $$w 1000000000000 -c 2>&1 >/dev/null); \
			echo "    10¹²: $$t s"; \
		else \
			echo "  Wheel $$w: 10⁹=$$med s ≥ 30s → SKIP 10¹²"; \
		fi; \
	done < /tmp/bench_times.$$; \
	rm -f /tmp/bench_times.$$

# PSVH reference generation
PSVH_DIR = psvh/primes

psvh-gen: $(TARGET)
	@mkdir -p $(PSVH_DIR)
	@for n in 3 4 5 6 7 8 9; do \
		target=$$((10**n)); \
		echo "Generating primes_1e$${n}.txt ..."; \
		./$(TARGET) --wheel 30 -o $(PSVH_DIR)/primes_1e$${n}.txt \
			--hash-output $(PSVH_DIR)/primes_1e$${n}.txt $$target >/dev/null 2>&1; \
		gzip -9 $(PSVH_DIR)/primes_1e$${n}.txt; \
	done
	@echo "Done. Files in $(PSVH_DIR)/"

psvh-verify: $(TARGET)
	@failed=0; \
	for f in $(PSVH_DIR)/primes_1e*.txt.gz; do \
		base=$$(basename "$$f" .gz); \
		dir=$$(dirname "$$f"); \
		hash=$$(cat "$$dir/$$base.sha256" | awk '{print $$1}'); \
		computed=$$(zcat "$$f" | ./$(TARGET) --verify-hash - 2>&1 | grep -o '[0-9a-f]\{64\}'); \
		if [ "$$hash" = "$$computed" ]; then \
			echo "  PASS: $$base"; \
		else \
			echo "  FAIL: $$base (expected $$hash, got $$computed)"; \
			failed=1; \
		fi; \
	done; \
	exit $$failed