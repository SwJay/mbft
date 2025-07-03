# bigquery.py
import argparse
import json
import os
import time
import logging
from decimal import Decimal
from google.cloud import bigquery
from google.cloud.exceptions import NotFound
from typing import List, Optional, Dict, Any, Tuple, Union # Add typing imports

# --- Configuration ---
# PROJECT_ID will be set from args or environment variable
DATASET_ID = "eth-query-456506.ethdat"
TABLE_ID = "eth-128"
# Use the same pattern as local.py for query file naming convention
# Add block_size (x) to the pattern
DEFAULT_QUERY_FILE_PATTERN = "query-x{block_size}-n{number}-w{window}-v{value_range:.1f}-{boolean_logic}-s{selectivity}-{agg_type}.json"
# Define the block number offset required by the user
BLOCK_NUMBER_OFFSET = 0

# --- Helper Functions (Copied/adapted from local.py) ---

def _validate_positive_int(value: str, arg_name: str) -> int:
    """Helper to validate positive integer arguments."""
    try:
        ival = int(value)
        if ival <= 0:
            raise argparse.ArgumentTypeError(f"{arg_name} must be a positive integer.")
        return ival
    except ValueError:
        raise argparse.ArgumentTypeError(f"{arg_name} must be an integer.")

def _validate_non_negative_int(value: str, arg_name: str) -> int:
    """Helper to validate non-negative integer arguments."""
    try:
        ival = int(value)
        if ival < 0:
            raise argparse.ArgumentTypeError(f"{arg_name} must be a non-negative integer.")
        return ival
    except ValueError:
        raise argparse.ArgumentTypeError(f"{arg_name} must be an integer.")

def _validate_range_float(value: str, arg_name: str) -> float:
    """Helper to validate float arguments within [0.0, 1.0]."""
    try:
        fval = float(value)
        if not (0.0 <= fval <= 1.0):
            raise argparse.ArgumentTypeError(f"{arg_name} must be between 0.0 and 1.0.")
        return fval
    except ValueError:
        raise argparse.ArgumentTypeError(f"{arg_name} must be a float.")

def _generate_query_filename(block_size: int, number: int, window: int, value_range: float, boolean_logic: str, selectivity: int, agg_type: str) -> str:
    """Generates the standard query filename based on parameters."""
    # Ensure boolean_logic and agg_type are formatted correctly before passing to format()
    # Note: boolean_logic and agg_type arrive already uppercased due to argparse type=str.upper
    # We need lowercase boolean and uppercase agg_type in the final filename.
    return DEFAULT_QUERY_FILE_PATTERN.format(
        block_size=block_size, # Add block_size here
        number=number,
        window=window,
        value_range=value_range,
        boolean_logic=boolean_logic.lower(), # Convert to lowercase HERE
        selectivity=selectivity,
        agg_type=agg_type # Already uppercase, pass directly
    )

def _parse_range(query_spec: Dict[str, Any], query_desc: str) -> Tuple[Optional[float], Optional[float]]:
    """Parses the 'range' field from a query specification."""
    low_amount, up_amount = None, None
    range_spec = query_spec.get("range")
    if range_spec:
        if isinstance(range_spec, list) and len(range_spec) > 0:
            # Expecting [[low, high]]
            if isinstance(range_spec[0], list) and len(range_spec[0]) == 2:
                val1, val2 = range_spec[0]
                try:
                    # Attempt conversion, allowing None if originally None
                    low_amount = float(val1) if val1 is not None else None
                    up_amount = float(val2) if val2 is not None else None
                    # Basic validation: low should not be greater than high if both are numbers
                    if low_amount is not None and up_amount is not None and low_amount > up_amount:
                         # Swap them
                         print(f"Warning: Query '{query_desc}' - range low value is greater than high value: {range_spec[0]}. Swapping them.")
                         low_amount, up_amount = up_amount, low_amount
                except (ValueError, TypeError):
                    # Reset on error
                    print(f"Warning: Query '{query_desc}' - Invalid range value (must be number or None): {range_spec[0]}. Skipping range filter.")
                    low_amount, up_amount = None, None
            else:
                 print(f"Warning: Query '{query_desc}' - Invalid 'range' inner list format (needs [low, high]): {range_spec[0]}. Skipping range filter.")
        else:
            print(f"Warning: Query '{query_desc}' - Invalid 'range' format (needs list): {range_spec}. Skipping range filter.")
    return low_amount, up_amount


def _parse_keyword_exp(query_spec: Dict[str, Any], query_desc: str) -> Tuple[str, List[str]]:
    """Parses the 'keyword_exp' field, returning logic operator and address list."""
    logic_operator = "OR"  # Default
    address_list_for_query: List[str] = []
    keyword_exp = query_spec.get("keyword_exp")

    if keyword_exp and isinstance(keyword_exp, dict):
        address_inputs: Optional[List[Any]] = None
        found_logic = False
        # Check for 'and' or 'or' keys (case-insensitive for robustness)
        logic_key = None
        if "and" in (key.lower() for key in keyword_exp):
            logic_key = next((k for k in keyword_exp if k.lower() == "and"), None)
            if logic_key and isinstance(keyword_exp[logic_key], list):
                 logic_operator, address_inputs, found_logic = "AND", keyword_exp[logic_key], True
        elif "or" in (key.lower() for key in keyword_exp):
            logic_key = next((k for k in keyword_exp if k.lower() == "or"), None)
            if logic_key and isinstance(keyword_exp[logic_key], list):
                 logic_operator, address_inputs, found_logic = "OR", keyword_exp[logic_key], True

        if not found_logic:
            print(f"Warning: Query '{query_desc}' - 'keyword_exp' must contain 'and' or 'or' key with a list value. Skipping keyword filter. {keyword_exp}")

        if found_logic and address_inputs is not None:
            extracted_addresses: List[str] = []
            for item in address_inputs:
                # Expecting {"input": "'address_value'"} or similar
                if isinstance(item, dict) and "input" in item and isinstance(item["input"], str):
                    # Strip quotes and whitespace
                    addr = item["input"].strip().strip("'\"")
                    # Basic check for valid address format
                    if addr and addr.startswith('0x'):
                        # Normalize to lower case
                        extracted_addresses.append(addr.lower())
                    else:
                        print(f"Warning: Query '{query_desc}' - Invalid or incorrectly formatted 'input' value in 'keyword_exp': {item}")
                else:
                    print(f"Warning: Query '{query_desc}' - Invalid item format in 'keyword_exp' (needs {{'input': 'addr'}}): {item}")

            if extracted_addresses:
                # Remove duplicates
                address_list_for_query = list(set(extracted_addresses))
            # Only warn if we expected addresses but got none
            elif found_logic:
                 print(f"Warning: Query '{query_desc}' - Failed to extract valid addresses from 'keyword_exp'.")

    return logic_operator, address_list_for_query

# --- Core BigQuery Functions ---

def execute_bigquery_query(
    client: bigquery.Client,
    query_spec: Dict[str, Any],
    query_desc: str,
    project_id: str # Pass project_id for table reference
) -> Tuple[Optional[Union[int, float]], float]:
    """
    Constructs and executes a single BigQuery query based on the spec.
    Amount range and aggregation now refer to the 'gas' column (gas limit).
    Returns (result, query_time_ms), where query_time_ms is the sum of compute times
    from the query plan stages, 0.0 if cached, or total job duration as fallback.
    """
    start_time = time.perf_counter()
    result_value: Optional[Union[int, float]] = None
    # Initialize BigQuery duration (used as fallback)
    bq_duration_ms: float = 0.0
    # This will hold the final reported time
    query_time_ms: float = 0.0
    # For error reporting
    final_sql = "N/A"
    # For error reporting
    params = []

    try:
        # 1. Extract parameters from query_spec
        agg_type = query_spec.get("agg_type", "").upper()
        begin_block = query_spec.get("start_blk")
        end_block = query_spec.get("end_blk")
        # low_amount/up_amount now refer to gas units based on user request
        low_amount, up_amount = _parse_range(query_spec, query_desc)
        logic_operator, address_list = _parse_keyword_exp(query_spec, query_desc)

        # Add COUNTDIST to valid types
        valid_agg_types = ["COUNT", "SUM", "MAX", "MIN", "COUNTDIST"]
        if agg_type not in valid_agg_types:
            # Return 0.0 for time on skipped query
            print(f"Warning: Query '{query_desc}' - Invalid aggregation type: {agg_type}. Skipping.")
            return None, 0.0

        # 2. Construct BigQuery SQL query with parameterization
        sql_parts = []
        # Reset params for this query
        params = []

        # SELECT clause based on agg_type, operating on the 'gas' column
        # The column to aggregate/filter on
        target_column = "gas"
        if agg_type == "COUNT":
            sql_parts.append("SELECT COUNT(*) as result")
        elif agg_type == "SUM":
            sql_parts.append(f"SELECT SUM({target_column}) as result")
        elif agg_type == "MAX":
             sql_parts.append(f"SELECT MAX({target_column}) as result")
        elif agg_type == "MIN":
             sql_parts.append(f"SELECT MIN({target_column}) as result")
        elif agg_type == "COUNTDIST":
            # Add COUNT(DISTINCT ...) for COUNTDIST
            sql_parts.append(f"SELECT COUNT(DISTINCT {target_column}) as result")

        # FROM clause - Use full path with project_id
        full_table_id = f"{project_id}.{DATASET_ID}.{TABLE_ID}"
        # Check if PROJECT_ID is already in DATASET_ID (it is for public datasets)
        if DATASET_ID.startswith(project_id):
             full_table_id = f"{DATASET_ID}.{TABLE_ID}"
        # If dataset ID is simple, prepend project
        elif '.' not in DATASET_ID:
             full_table_id = f"{project_id}.{DATASET_ID}.{TABLE_ID}"
        # Assume DATASET_ID includes project if it contains '.'
        else:
             full_table_id = f"{DATASET_ID}.{TABLE_ID}"

        sql_parts.append(f"FROM `{full_table_id}`")

        # WHERE clause
        where_clauses = []
        # Block filter
        if begin_block is not None:
            try:
                # Apply the offset before creating the parameter
                adjusted_begin_block = int(begin_block) + BLOCK_NUMBER_OFFSET
                where_clauses.append("block_number >= @start_blk")
                params.append(bigquery.ScalarQueryParameter("start_blk", "INT64", adjusted_begin_block))
            except (ValueError, TypeError):
                 print(f"Warning: Query '{query_desc}' - Invalid start_blk value {begin_block}. Ignoring lower block bound filter.")
        if end_block is not None:
            try:
                # Apply the offset before creating the parameter
                adjusted_end_block = int(end_block) + BLOCK_NUMBER_OFFSET
                where_clauses.append("block_number <= @end_blk")
                params.append(bigquery.ScalarQueryParameter("end_blk", "INT64", adjusted_end_block))
            except (ValueError, TypeError):
                print(f"Warning: Query '{query_desc}' - Invalid end_blk value {end_block}. Ignoring upper block bound filter.")

        # Amount (Gas Limit) filter
        if low_amount is not None:
            # Compare directly against the target column ('gas')
            where_clauses.append(f"{target_column} >= @low_amount")
            # Use FLOAT64 for parameter type to match _parse_range output and allow flexible ranges
            params.append(bigquery.ScalarQueryParameter("low_amount", "FLOAT64", low_amount))
        if up_amount is not None:
             # Compare directly against the target column ('gas')
             where_clauses.append(f"{target_column} <= @up_amount")
             # Use FLOAT64 for parameter type
             params.append(bigquery.ScalarQueryParameter("up_amount", "FLOAT64", up_amount))

        # Address filter
        if address_list:
            address_param_name = "address_list"
            params.append(bigquery.ArrayQueryParameter(address_param_name, "STRING", address_list))
            # Use OR logic: transaction involves any of the addresses in from_address OR to_address
            address_clause = f"(from_address IN UNNEST(@{address_param_name}) OR to_address IN UNNEST(@{address_param_name}))"
            if logic_operator == "AND":
                 # Note: This interpretation of AND (both from AND to must be in the list for a single tx)
                 # is different from MongoDB's $all logic on an array field.
                 # We proceed with the OR logic as it's a closer match to the likely intent.
                 # If a strict AND (both from AND to in list) were needed, the clause would be:
                 # address_clause = f"(from_address IN UNNEST(@{address_param_name}) AND to_address IN UNNEST(@{address_param_name}))"
                 print(f"Info: Query '{query_desc}' - 'AND' logic has ambiguity with BigQuery from/to address structure. Current implementation uses 'OR' logic (either address matches).")
            where_clauses.append(address_clause)


        if where_clauses:
            sql_parts.append("WHERE " + "\n  AND ".join(where_clauses))

        final_sql = "\n".join(sql_parts)
        # Optional: print SQL for debugging
        # print(f"  Executing SQL:\n{final_sql}")
        # print(f"  Parameters: {[(p.name, p.value) for p in params]}")

        # 3. Configure and run the query job
        job_config = bigquery.QueryJobConfig(query_parameters=params)
        query_job = client.query(final_sql, job_config=job_config)

        # 4. Wait for results and extract
        # .result() waits for job completion.
        rows_iterator = query_job.result()

        # --- Attempt to reload the job object to ensure statistics are populated ---
        try:
            # Reload the job object using its ID and location
            # print(f"  Debug Info: Reloading job object for latest status (Job ID: {query_job.job_id})...")
            query_job = client.get_job(query_job.job_id, location=query_job.location)
            # print(f"  Debug Info: Job object reloaded.")
        except Exception as reload_e:
             print(f"  Debug Info: Error reloading job object: {reload_e} - Will proceed with original job object.")
        # ---------------------------------------------------------------------------

        # --- Get BigQuery execution time (total job duration fallback) ---
        try:
            if query_job.started and query_job.ended:
                bq_duration = query_job.ended - query_job.started
                bq_duration_ms = bq_duration.total_seconds() * 1000
            else:
                 # This might happen if the job failed very early or info is unavailable
                 # Set fallback to 0 if timestamps unavailable
                 print("  Query Info: BigQuery job start/end timestamps unavailable.")
                 bq_duration_ms = 0.0
        except Exception as time_e:
             # Set fallback to 0 on error
             print(f"  Query Info: Error getting BigQuery job timestamps: {time_e}")
             bq_duration_ms = 0.0
        # ------------------------------------

        # Extract the single result value
        # BQ returns an iterator even for aggregate queries that produce one row
        for row in rows_iterator:
            # Access result by alias 'result'
            result_value = row.result
            # Stop after the first (and only) row
            break

        # Handle cases where query returns no rows or NULL result
        # This can happen for MIN/MAX/SUM on an empty filtered set. COUNT returns 0.
        if result_value is None:
            # Sum of empty set is 0
            if agg_type in ["COUNT", "SUM", "COUNTDIST"]:
                 result_value = 0
             # For MIN/MAX, None is the appropriate result for an empty set

        # Convert BigQuery's Decimal type to float for consistency
        if isinstance(result_value, Decimal):
             result_value = float(result_value)

        # --- Calculate Query Time (Compute or Fallback) ---
        total_compute_ms = 0.0
        plan_available = False
        try:
            # Check cache hit first
            if getattr(query_job, 'cache_hit', False):
                # Cached query has 0 compute time
                print("  Query Info: Query result from cache, compute time recorded as 0.0 ms.")
                query_time_ms = 0.0
                # Mark plan as "handled" (cached)
                plan_available = True
            # If not cache hit, check for query_plan directly on the job object
            elif hasattr(query_job, 'query_plan') and query_job.query_plan:
                plan_available = True
                print("  --- Query Stage Compute Time Details ---")
                # Access query_plan directly
                for stage in query_job.query_plan:
                    # Safely get timing attributes, defaulting to 0 if None or missing
                    # wait_ms = getattr(stage, 'wait_ms_avg', 0) or 0
                    # read_ms = getattr(stage, 'read_ms_avg', 0) or 0
                    compute_ms = getattr(stage, 'compute_ms_avg', 0) or 0
                    # write_ms = getattr(stage, 'write_ms_avg', 0) or 0
                    stage_name = getattr(stage, 'name', 'Unknown Stage')
                    total_compute_ms += compute_ms
                    # Simplified print focusing on compute time
                    print(f"    Stage: {stage_name:<15} | Compute: {compute_ms:>8.2f} ms")
                query_time_ms = total_compute_ms
                print(f"  --- Total Compute Time: {total_compute_ms:.2f} ms ---")
            # Fallback if not cache hit and no query_plan found (should be unusual)
            else:
                 # Use total duration as fallback
                 print(f"  Query Info: Query plan not found. Using total job duration {bq_duration_ms:.2f} ms as fallback.")
                 query_time_ms = bq_duration_ms
        except Exception as plan_e:
            # Use total duration as fallback
            print(f"  Query Info: Error extracting query plan timing: {plan_e}. Using total job duration {bq_duration_ms:.2f} ms as fallback.")
            query_time_ms = bq_duration_ms
        # ----------------------------------------------------------

        # Print cost/bytes processed info if available
        # Access these directly from query_job as well
        try:
             # Use getattr for safety, though these should exist on a completed job
             bytes_billed = getattr(query_job, 'total_bytes_billed', None)
             bytes_processed = getattr(query_job, 'total_bytes_processed', None)

             if bytes_billed is not None:
                  gb_billed = bytes_billed / (1024 ** 3)
                  # Note: Cost depends on region and pricing tier. This is an estimate based on bytes.
                  # Using a rough estimate of $6.25 USD per TB (as of late 2023/early 2024 for on-demand)
                  cost_usd = (bytes_billed / (1024**4)) * 6.25
                  print(f"  Query Info: Processed {gb_billed:.4f} GB, Estimated cost: ${cost_usd:.4f} USD")
             elif bytes_processed is not None:
                  gb_processed = bytes_processed / (1024 ** 3)
                  print(f"  Query Info: Processed {gb_processed:.4f} GB (Read from cache or materialized view, cost may be $0)")

        except Exception as billing_e:
            print(f"  Query Info: Could not retrieve billing information: {billing_e}")


    except NotFound as e:
        # Ensure result is None on error
        print(f"Error: Query '{query_desc}' failed - BigQuery table or dataset not found: {e}")
        result_value = None
        # Assign 0 time on fatal error
        query_time_ms = 0.0
    except Exception as e:
        print(f"Error: Unexpected error executing query '{query_desc}': {e}")
        print(f"  Failed SQL:\n{final_sql}")
        # Safely print parameters, handling different types
        param_repr = []
        for p in params:
            if isinstance(p, bigquery.ScalarQueryParameter):
                param_repr.append((p.name, p.value))
            # Use array_values
            elif isinstance(p, bigquery.ArrayQueryParameter):
                param_repr.append((p.name, p.array_values))
            # Fallback
            else:
                 param_repr.append((p.name, f"Unknown parameter type: {type(p)}"))
        print(f"  Failed Parameters: {param_repr}")
        # Ensure result is None on error
        result_value = None
        # Assign 0 time on fatal error
        query_time_ms = 0.0

    # Return the calculated or fallback query time
    return result_value, query_time_ms

def execute_queries_from_file_bq(query_file: str, client: bigquery.Client, project_id: str, target_table_id: Optional[str] = None, placeholder_table_id: Optional[str] = None) -> Optional[float]:
    """ Reads query definitions from a JSON file and executes them using BigQuery.
        Returns the average BigQuery compute time (sum of stages' compute_ms) in ms
        for successfully executed queries, or None if none were executed.
        Returns 0.0 for cached queries or if query plan is unavailable (uses total duration as fallback in that case).
    """
    try:
        if not os.path.exists(query_file):
            # Return 0.0 if file not found
            print(f"Error: Query file '{query_file}' not found.")
            return 0.0
        with open(query_file, 'r') as f:
            queries_data = json.load(f)
    except json.JSONDecodeError as e:
        # Return 0.0 on parse error
        print(f"Error: Failed to parse query file '{query_file}': {e}")
        return 0.0
    except IOError as e:
        # Return 0.0 on read error
        print(f"Error: Error reading query file '{query_file}': {e}")
        return 0.0
    # Catch other potential errors during file read/load
    except Exception as e:
        # Return 0.0 on other error
        print(f"An unknown error occurred while reading or parsing the query file: {e}")
        return 0.0

    print(f"\n--- Executing BigQuery queries from '{query_file}' ({project_id}) ---")
    if not isinstance(queries_data, list):
        # Return 0.0 if format is wrong
        print(f"Error: The top-level structure of query file '{query_file}' must be a JSON array/list.")
        return 0.0

    # Update variable names for clarity
    total_query_compute_time_ms = 0.0
    executed_query_count = 0
    # Count queries that ran without error AND returned a value
    successful_query_count = 0

    for i, query_spec in enumerate(queries_data):
        query_num = i + 1
        print(f"\n--- Query {query_num} ---")
        if not isinstance(query_spec, dict):
            print(f"Warning: Skipping invalid query entry {query_num} (not a dictionary): {query_spec}")
            continue

        # Count that we attempted to execute it
        executed_query_count += 1
        query_desc = query_spec.get("query_desc", f"Unnamed Query {query_num}")

        # --- Print query details (similar to local.py) ---
        agg_type = query_spec.get("agg_type", "N/A").upper()
        begin_block = query_spec.get("start_blk")
        end_block = query_spec.get("end_blk")
        # Values from range now refer to Gas Limit units
        low_gas_limit, up_gas_limit = _parse_range(query_spec, query_desc)
        logic_operator, address_list_for_query = _parse_keyword_exp(query_spec, query_desc)

        print(f"Description: {query_desc}")
        print(f"  Aggregation Type: {agg_type}")
        print(f"  Block Range: {begin_block if begin_block is not None else 'N/A'} - {end_block if end_block is not None else 'N/A'}")
        # Update print statement to reflect gas limit
        print(f"  Gas Limit Range: {low_gas_limit if low_gas_limit is not None else 'N/A'} - {up_gas_limit if up_gas_limit is not None else 'N/A'} (units)")
        # Info message about AND/OR handled in execute function
        print(f"  Address Logic: {logic_operator}")
        print(f"  Address Filter: {address_list_for_query if address_list_for_query else 'None'}")


        # --- Execute the BigQuery query ---
        # Get the specific query time (compute, cache=0, or fallback=duration)
        result, single_query_time_ms = execute_bigquery_query(client, query_spec, query_desc, project_id)

        # --- Accumulate stats ---
        # Accumulate the returned time (compute, 0 for cache, or duration fallback)
        total_query_compute_time_ms += single_query_time_ms
        # Count only successfully executed queries that returned a non-error value
        if result is not None:
             successful_query_count += 1

        # --- Print result ---
        # Format float results nicely (SUM/MAX/MIN of gas might be large ints or potentially floats if BQ processes them as such)
        if isinstance(result, float):
            # Check if it's effectively an integer
            if result.is_integer():
                print(f"Result: {int(result)} (Gas Units)")
            else:
                # Keep float formatting if it has decimals, though unlikely for gas
                print(f"Result: {result:.2f} (Gas Units)")
        elif isinstance(result, int):
             print(f"Result: {result} (Gas Units)")
        else:
            # Handle None result explicitly (e.g., for MIN/MAX on empty set or error)
            print(f"Result: {result if result is not None else 'No matching documents or query failed'}")
        # Update print label for time
        print(f"BigQuery Query Time (Compute/Cache/Fallback): {single_query_time_ms:.2f} ms")


    # --- Print overall stats ---
    print("\n--- BigQuery JSON Query Execution Finished ---")
    if executed_query_count > 0:
        # Calculate average based on the accumulated time
        average_query_time_ms = total_query_compute_time_ms / executed_query_count
        print(f"Attempted to execute {executed_query_count} queries.")
        print(f"Successfully returned results for {successful_query_count} queries.")
        # Update print labels for time
        print(f"Total BigQuery Query Time (Compute/Cache/Fallback): {total_query_compute_time_ms:.2f} ms.")
        print(f"Average BigQuery Query Time per query (Compute/Cache/Fallback): {average_query_time_ms:.2f} ms.")
        # Return the calculated average time
        return average_query_time_ms
    else:
        print("Query file was empty or all entries were invalid, no queries were executed.")
        # Return 0.0 if no queries were executed
        return 0.0


# --- Main execution block ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Query bigquery-public-data.crypto_ethereum.transactions data using Google BigQuery.\n"
                    "Requires a valid GCP Project ID and authentication (e.g., 'gcloud auth application-default login').",
        formatter_class=argparse.RawTextHelpFormatter
    )

    # Although only 'query' is implemented, keep structure for potential future actions
    # parser.add_argument('--action', default='query', choices=['query'], help='Action to perform (currently only query is supported)')

    # parser.add_argument('--project-id', type=str, default=os.environ.get("GCP_PROJECT"),
    #                      help='Your Google Cloud Project ID (or set GCP_PROJECT environment variable). This is required.')
    parser.add_argument('--project-id', type=str, default="eth-query-456506", help='Your Google Cloud Project ID. This is required.')

    # Arguments for locating the query file (similar to local.py's query action)
    # Make group optional
    query_group = parser.add_mutually_exclusive_group(required=False)
    query_group.add_argument('-f', '--file', type=str, default="query/query_1.json",
                             help='Specify the path to the query JSON file to execute (takes precedence over auto-generated filename).')

    # Add arguments for filename generation (if --file is not provided)
    # Reuse validation helpers and defaults from local.py
    parser.add_argument('-n', '--number', type=lambda v: _validate_positive_int(v, '-n/--number'), default=100,
                              help="Number of queries (used for filename generation, default: 100)")
    parser.add_argument('-w', '--window', type=lambda v: _validate_non_negative_int(v, '-w/--window'), default=800,
                              help="Time window size (used for filename generation, default: 800)")
    # Updated help text
    parser.add_argument('-v', '--value-range', type=lambda v: _validate_range_float(v, '-v/--value-range'), default=0.5,
                              help="Gas Limit range ratio [0.0-1.0] (used for filename generation, default: 0.5)")
    parser.add_argument('-b', '--boolean', choices=['AND', 'OR'], default='OR', type=str.upper,
                              help="Boolean logic (used for filename generation, default: OR)")
    parser.add_argument('-s', '--selectivity', type=lambda v: _validate_positive_int(v, '-s/--selectivity'), default=2,
                              help="Keyword selectivity (used for filename generation, default: 2)")
    # Updated help text & choices
    # Updated choices to include MIN
    parser.add_argument('-a', '--agg-type', choices=['MAX', 'MIN', 'COUNT', 'SUM', 'COUNTDIST'], default='COUNT', type=str.upper,
                              help="Aggregation type (targets Gas Limit, default: COUNT)")

    # Add --case argument mirroring local.py
    parser.add_argument('-c', '--case', type=lambda v: _validate_non_negative_int(v, '-c/--case'), default=1,
                            help="Select execution mode:\n" \
                                "  0: Single query file (specified by --file or generated from params)\n" \
                                "  1: Vary time window (window), test both AND and OR logic\n" \
                                "  2: Vary Value Range (value-range)\n" \
                                "  3: Vary keyword number (selectivity)\n" \
                                "  4: Vary Block Size (affects filename only)\n" \
                                "  5: Vary aggregation type (agg-type)")

    try:
        args = parser.parse_args()
    # Exit gracefully if validation helpers raise ArgumentTypeError
    except argparse.ArgumentTypeError as e:
         parser.error(str(e))

    # --- Validate Project ID ---
    if not args.project_id:
         parser.error("Error: Google Cloud Project ID is required.\nPlease use the --project-id argument or set the GCP_PROJECT environment variable.")
    current_project_id = args.project_id


    # --- Initialize BigQuery Client ---
    bq_client: Optional[bigquery.Client] = None
    try:
        print(f"Initializing BigQuery client (Project: {current_project_id})...")
        # Enable debug logging for google-cloud libs
        logging.basicConfig(level=logging.DEBUG)
        # Add custom log
        print("DEBUG: Preparing to call bigquery.Client()")
        bq_client = bigquery.Client(project=current_project_id)
        # Add custom log
        print("DEBUG: bigquery.Client() call completed")
        print("Testing BigQuery connection...")
        # Add custom log
        print("DEBUG: Preparing to call bq_client.list_datasets()")
        # Simple call to test
        list(bq_client.list_datasets(max_results=1))
        # Add custom log
        print("DEBUG: bq_client.list_datasets() call completed")
        print("BigQuery client initialized and connected successfully.")
        # Optional: Reset logging level if too verbose later
        logging.basicConfig(level=logging.INFO)
    except Exception as e:
        print(f"\nError: Failed to initialize or connect to BigQuery client: {e}")
        print("\nPlease check the following:")
        print(f"  1. Is the GCP Project ID '{current_project_id}' correct?")
        print("  2. Have you authenticated using 'gcloud auth application-default login'?")
        print(f"  3. Is the BigQuery API enabled in project '{current_project_id}'?")
        print("  4. Do your authenticated credentials have sufficient permissions to access BigQuery resources?")
        # Ensure client is None if init fails
        bq_client = None
        # exit(1) # Optionally exit, or let the subsequent checks handle None client

    # --- Execute Queries based on --case ---
    # Proceed only if the client was initialized successfully
    if bq_client:
        print(f"\n--- Selected Action: Execute Queries (Case: {args.case}) ---")

        # Default parameters (similar to local.py, add block_size)
        # Example default block size
        default_block_size = 128
        # Use provided arg or its default
        default_number = args.number
        default_window = args.window
        default_value_range = args.value_range
        default_boolean = args.boolean
        default_selectivity = args.selectivity
        default_agg_type = args.agg_type

        # Define the base output directory for BigQuery results
        bq_res_dir = "res_bq"
        # Assuming query files are in subdirs like 'query/window_and/', 'query/value_range/' etc.
        bq_query_dir_base = "query_bq"

        if args.case == 0:
            print("Mode: Single query file")
            query_file_to_use: Optional[str] = None
            # User specified a file directly
            if args.file:
                query_file_to_use = args.file
                print(f"Using user-specified query file: '{query_file_to_use}'")
                if not os.path.exists(query_file_to_use):
                     # Prevent execution
                     print(f"Error: Specified query file '{query_file_to_use}' does not exist.")
                     query_file_to_use = None
            # --file not provided, generate filename from other args
            else:
                 print("Parameter --file not specified, generating filename based on other parameters...")
                 # Generate filename using defaults or provided args
                 query_file_to_use = _generate_query_filename(
                     # Need a block size for filename
                     default_block_size,
                     args.number, args.window, args.value_range,
                     args.boolean, args.selectivity, args.agg_type
                 )
                 print(f"Auto-generated query filename: '{query_file_to_use}' (assuming in . or {bq_query_dir_base}/ directory)")
                 # Check existence (might need adjustment based on where generated files are expected)
                 if not os.path.exists(query_file_to_use) and not os.path.exists(os.path.join(bq_query_dir_base, query_file_to_use)):
                      # Prevent execution
                      print(f"Warning: Auto-generated query file '{query_file_to_use}' not found in the current directory or '{bq_query_dir_base}/'.")
                      print("Please ensure the file exists.")
                      query_file_to_use = None
                 elif os.path.exists(os.path.join(bq_query_dir_base, query_file_to_use)) and not os.path.exists(query_file_to_use):
                     # Use path within query dir
                     query_file_to_use = os.path.join(bq_query_dir_base, query_file_to_use)

            if query_file_to_use:
                execute_queries_from_file_bq(query_file_to_use, bq_client, current_project_id, target_table_id=None, placeholder_table_id=None)
            else:
                 print("--- No valid query file found, no queries executed ---")

        elif args.case == 1:
            print("Mode: Vary Time Window (W), AND/OR Boolean Logic (B)")
            booleans = ['AND', 'OR']
            windows = [400, 800, 1200, 1600, 2000]

            for b in booleans:
                logic_lower = b.lower()
                # Define output file path within res_bq directory
                output_file = os.path.join(bq_res_dir, f'window_{logic_lower}.csv')
                output_dir = os.path.dirname(output_file)
                os.makedirs(output_dir, exist_ok=True)
                print(f"Results will be written to: {output_file}")

                # Define query file directory for this sub-case
                query_subdir = f"window_{logic_lower}"
                query_dir = os.path.join(bq_query_dir_base, query_subdir)

                with open(output_file, 'w') as f:
                    # Update CSV header
                    f.write("time window,\tavg compute time (ms)\n")

                    for w in windows:
                        query_file_name = _generate_query_filename(
                            default_block_size, default_number, w, default_value_range, b, default_selectivity, default_agg_type
                        )
                        query_file_path = os.path.join(query_dir, query_file_name)
                        # Default to N/A
                        avg_time_ms: Union[float, str] = "N/A"
                        if os.path.exists(query_file_path):
                            print(f"  Executing: {query_file_path}")
                            query_time_result = execute_queries_from_file_bq(query_file_path, bq_client, current_project_id, target_table_id=None, placeholder_table_id=None)
                            # Should return 0.0 if no queries ran
                            if query_time_result is not None:
                                avg_time_ms = f"{query_time_result:.2f}"
                        else:
                            print(f"  Warning: Query file not found, skipping: {query_file_path}")

                        f.write(f"{w}\t{avg_time_ms}\n")

        elif args.case == 2:
            print("Mode: Vary Value Range (V), Selectivity=0")
            value_ranges = [0.1, 0.3, 0.5, 0.7, 0.9]

            output_file = os.path.join(bq_res_dir, 'value_range.csv')
            output_dir = os.path.dirname(output_file)
            os.makedirs(output_dir, exist_ok=True)
            print(f"Results will be written to: {output_file}")

            query_subdir = "value_range"
            query_dir = os.path.join(bq_query_dir_base, query_subdir)

            with open(output_file, 'w') as f:
                # Update CSV header
                f.write("value range,\tavg compute time (ms)\n")

                for v in value_ranges:
                    query_file_name = _generate_query_filename(
                        default_block_size, default_number, default_window, v, default_boolean, default_selectivity, default_agg_type
                    )
                    query_file_path = os.path.join(query_dir, query_file_name)
                    avg_time_ms: Union[float, str] = "N/A"
                    if os.path.exists(query_file_path):
                        print(f"  Executing: {query_file_path}")
                        query_time_result = execute_queries_from_file_bq(query_file_path, bq_client, current_project_id, target_table_id=None, placeholder_table_id=None)
                        if query_time_result is not None:
                             avg_time_ms = f"{query_time_result:.2f}"
                    else:
                        print(f"  Warning: Query file not found, skipping: {query_file_path}")

                    # Format float value
                    f.write(f"{v:.1f}\t{avg_time_ms}\n")

        elif args.case == 3:
            # Represents selectivity
            print("Mode: Vary Keyword Number (S)")
            keywords = [2, 4, 8, 16, 32]

            output_file = os.path.join(bq_res_dir, 'keyword.csv')
            output_dir = os.path.dirname(output_file)
            os.makedirs(output_dir, exist_ok=True)
            print(f"Results will be written to: {output_file}")

            query_subdir = "keyword"
            query_dir = os.path.join(bq_query_dir_base, query_subdir)

            with open(output_file, 'w') as f:
                # Update CSV header
                f.write("keyword number (selectivity),\tavg compute time (ms)\n")

                for k in keywords:
                    query_file_name = _generate_query_filename(
                        default_block_size, default_number, default_window, default_value_range, default_boolean, k, default_agg_type
                    )
                    query_file_path = os.path.join(query_dir, query_file_name)
                    avg_time_ms: Union[float, str] = "N/A"
                    if os.path.exists(query_file_path):
                         print(f"  Executing: {query_file_path}")
                         query_time_result = execute_queries_from_file_bq(query_file_path, bq_client, current_project_id, target_table_id=None, placeholder_table_id=None)
                         if query_time_result is not None:
                              avg_time_ms = f"{query_time_result:.2f}"
                    else:
                        print(f"  Warning: Query file not found, skipping: {query_file_path}")

                    f.write(f"{k}\t{avg_time_ms}\n")

        elif args.case == 4:
            # Note: Block size only affects the *naming* of the query file for BigQuery,
            # as the actual data source is the fixed BigQuery table.
            print("Mode: Vary Block Size (X) (affects query filename only)")
            block_sizes = [32, 64, 128, 256, 512]

            output_file = os.path.join(bq_res_dir, 'blk_size.csv')
            output_dir = os.path.dirname(output_file)
            os.makedirs(output_dir, exist_ok=True)
            print(f"Results will be written to: {output_file}")

            query_subdir = "blk_size"
            query_dir = os.path.join(bq_query_dir_base, query_subdir)

            with open(output_file, 'w') as f:
                # Update CSV header
                f.write("block size (filename),\tavg compute time (ms)\n")

                for b_size in block_sizes:
                    query_file_name = _generate_query_filename(
                        b_size, default_number, default_window, default_value_range, default_boolean, default_selectivity, default_agg_type
                    )
                    query_file_path = os.path.join(query_dir, query_file_name)

                    # --- Modification Start ---
                    # Construct the target table name and ID based on block size
                    target_table_name = f"eth-{b_size}"
                    # ASSUMPTION: Dataset ID is 'eth_data'. Replace if different.
                    # ASSUMPTION: current_project_id holds your GCP project ID.
                    target_table_id = f"{current_project_id}.eth_data.{target_table_name}"
                    # Define the placeholder used in the SQL query files
                    # ASSUMPTION: Placeholder is 'eth-placeholder' in the 'eth_data' dataset.
                    placeholder_table_id = f"{current_project_id}.eth_data.eth-placeholder"
                    # --- Modification End ---

                    avg_time_ms: Union[float, str] = "N/A"
                    if os.path.exists(query_file_path):
                         print(f"  Executing query: {query_file_path} targeting table: {target_table_id}")
                         # --- Modification Start ---
                         # Pass target_table_id and placeholder_table_id to the execution function.
                         # NOTE: execute_queries_from_file_bq needs to be updated to accept
                         #       these arguments and perform the replacement internally.
                         query_time_result = execute_queries_from_file_bq(
                             query_file_path,
                             bq_client,
                             current_project_id,
                             # Pass placeholder too
                             target_table_id=target_table_id,
                             placeholder_table_id=placeholder_table_id
                         )
                         # --- Modification End ---
                         if query_time_result is not None:
                              avg_time_ms = f"{query_time_result:.2f}"
                    else:
                        print(f"  Warning: Query file not found, skipping: {query_file_path}")

                    f.write(f"{b_size}\t{avg_time_ms}\n")

        elif args.case == 5:
            print("Mode: Vary Aggregation Type (A)")
            # Use the aggregation types supported by the BigQuery implementation
            agg_types = ['MAX', 'COUNT', 'SUM', 'COUNTDIST']

            output_file = os.path.join(bq_res_dir, 'agg_type.csv')
            output_dir = os.path.dirname(output_file)
            os.makedirs(output_dir, exist_ok=True)
            print(f"Results will be written to: {output_file}")

            query_subdir = "agg_type"
            query_dir = os.path.join(bq_query_dir_base, query_subdir)

            with open(output_file, 'w') as f:
                # Update CSV header
                f.write("agg type,\tavg compute time (ms)\n")

                for a_type in agg_types:
                    query_file_name = _generate_query_filename(
                        default_block_size, default_number, default_window, default_value_range, default_boolean, default_selectivity, a_type
                    )
                    query_file_path = os.path.join(query_dir, query_file_name)
                    avg_time_ms: Union[float, str] = "N/A"
                    if os.path.exists(query_file_path):
                         print(f"  Executing: {query_file_path}")
                         query_time_result = execute_queries_from_file_bq(query_file_path, bq_client, current_project_id, target_table_id=None, placeholder_table_id=None)
                         if query_time_result is not None:
                              avg_time_ms = f"{query_time_result:.2f}"
                    else:
                        print(f"  Warning: Query file not found, skipping: {query_file_path}")

                    f.write(f"{a_type}\t{avg_time_ms}\n")

        else:
             print(f"Error: Unknown case value: {args.case}")

    else:
         print("\n--- BigQuery client initialization failed, cannot execute any queries ---")


    print("\nScript execution finished.")
