# README


## Usage

```
positional arguments:
    {load,query,delete,create-query} choose the operation to perform (load, query, create-query, delete)

    load                load data eth-{block_size}.dat to databse according to block size
        
        -x BLOCK_SIZE, --block-size BLOCK_SIZE
                                block size of data file (e.g., 128 load eth-128.dat)
    
    query               query from json files

        -f FILE, --file FILE    name for the query file

        -c CASE, --case CASE    case 0: single query following parameters
                                case 1: vary time window, both OR & AND
                                case 2: vary value range
                                case 3: vary keyword number
                                case 4: vary block size
                                case 5: vary aggregation type

    delete              delete entire db

    create-query        create query json file according to parameters

        -x BLOCK_SIZE, --block-size BLOCK_SIZE
                                database with differnt block size to create queries from

        -n NUMBER, --number NUMBER
                                number of queries

        -w WINDOW, --window WINDOW
                                time window size

        -v VALUE_RANGE, --value-range VALUE_RANGE
                                numerical value range

        -b {AND,OR}, --boolean  {AND,OR}
                                boolean function for keywords

        -s SELECTIVITY, --selectivity SELECTIVITY
                                number of keywords

        -a {MAX,COUNT,SUM,COUNTDIST}, --agg-type {MAX,COUNT,SUM,COUNTDIST}
                                aggregate type

        -c CASE, --case CASE    case 0: single query following parameters
                                case 1: vary time window, both OR & AND
                                case 2: vary value range
                                case 3: vary keyword number
                                case 4: vary block size
                                case 5: vary aggregation type
```

Example of usage:

```
python mongodb.py load -x 128   # create db
python mongodb.py create-query -x 128 -n 100 -w 800 -v 0.5 -b 'OR' -s 2 -a 1    # create customed query
python mongodb.py query -f path/to/query_file.json
```