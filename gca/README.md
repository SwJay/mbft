# README


## Usage

Requires GMP and PBC library.

```
./mbft [-option]

optional arguments:
    -s, --seed              Seed for randomness in Hash & BF
    -q, --universe-bits     Universe bits
    -v, --vmax-bits         Max value
    -i, --maxid-bits        Max id
    -d, --txn-data-path     Path to transaction data
    -p, --pbc-param-path    Path to PBC parameters
    -k, --key-dir           Path to the public key
    -a, --agg-type          Query aggregation type: COUNT, SUM, MAX, MIN, AVG
    -l, --lower             Lower bound of amount
    -u, --upper             Upper bound of amount
    -b, --begin             Start block
    -e, --end               End block
    -r, --repeats           Number of repeats
    -c, --test-case         Test cases
                            * 0. GenKeys:     Generate keys for the chain
                            * 1. Single shot:  Build a chain, do single query & verification
                            * 2. Vary window: test performance under verying time window
                            * 3. Vary value range: test performance under varying value range
                            * 4. Vary block size: test performance under varying block size
                            * 5. Vary aggregate type: test performance under varying aggregate type
```

Example of usage:

```
./gca_tree -c 0 -q {universe_bits} -k path/to/key -p path/to/param
./gca_tree -q 14 -v 9 -i 5 -c 1 -d path/to/dataset -k path/to/key -a 0 -l 1000 -u 2000 -b 100 -e 200    # run single query
./gca_tree -q 14 -v 9 -i 5 -c 2 -d path/to/dataset -k path/to/key   # Other test cases
```
