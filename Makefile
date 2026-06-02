CFLAGS = -O3 -pipe -Wall -Wextra
LDFLAGS =
TARGET = fastsieve

.PHONY: all clean test

all: $(TARGET)

$(TARGET): fastsieve.c
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^

clean:
	rm -f $(TARGET)

test: $(TARGET)
	@got=$$(./$(TARGET) 1000 -c 2>&1 | sed -n 's/Found \([0-9]*\) primes.*/\1/p'); \
	[ "$$got" = "168" ] || { echo "FAIL: $(TARGET) 1000 -c => $$got, expected 168"; exit 1; }; \
	echo "PASS: π(1000) = 168"
	@got=$$(./$(TARGET) 1000000 -c 2>&1 | sed -n 's/Found \([0-9]*\) primes.*/\1/p'); \
	[ "$$got" = "78498" ] || { echo "FAIL: $(TARGET) 1000000 -c => $$got, expected 78498"; exit 1; }; \
	echo "PASS: π(10⁶) = 78498"
	@got=$$(./$(TARGET) 1000000000 -c 2>&1 | sed -n 's/Found \([0-9]*\) primes.*/\1/p'); \
	[ "$$got" = "50847534" ] || { echo "FAIL: $(TARGET) 1000000000 -c => $$got, expected 50847534"; exit 1; }; \
	echo "PASS: π(10⁹) = 50847534"
