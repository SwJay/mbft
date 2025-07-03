# README

dataset id should be modified according to your google bigquery table

## Usage

```
    --project-id PROJECT_ID
                            Google Cloud project ID.

    -f FILE, --file FILE  query file

    -n NUMBER, --number NUMBER
                            number of quries

    -w WINDOW, --window WINDOW
                            time window

    -v VALUE_RANGE, --value-range VALUE_RANGE
                            amount range

    -b {AND,OR}, --boolean {AND,OR}
                            boolean function

    -s SELECTIVITY, --selectivity SELECTIVITY
                            number of keywords

    -a {MAX,MIN,COUNT,SUM,COUNTDIST}, --agg-type {MAX,MIN,COUNT,SUM,COUNTDIST}
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
python bigquery.py -f   # query selected file
```