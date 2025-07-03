# local.py
import pymongo
import re
import os
import argparse
import json
import random
import time
from decimal import Decimal, InvalidOperation
from typing import List, Optional, Dict, Any, Tuple, Union # Add typing imports
from pymongo.collection import Collection # Import Collection for type hinting

# --- Helper Function for Collection Name ---
def get_collection_name(data_file_path: str) -> str:
    """Generates a MongoDB collection name from a data file path."""
    base_name = os.path.basename(data_file_path)
    name_part, _ = os.path.splitext(base_name)
    # Basic sanitization: replace characters potentially unsafe for collection names or inconsistent style
    # Replace ., $, -, \\, /, space, ' with _
    safe_name_part = re.sub(r'[.$\-\\\\/ \']', '_', name_part)
    return f"transactions_{safe_name_part}"

# --- Configuration ---
MONGO_URI = "mongodb://localhost:27017/"
DB_NAME = "blockchain_db"
DATA_FILE = "../data/eth.dat" # Default data file, will be overridden by args for load/create-query

# Regular expression to parse each line of the data file
LINE_REGEX = re.compile(r"^\s*(\d+)\s+\[\s*([\d.]+)\s*\]\s+\{\s*(.*?)\s*\}\s*$")

# --- Helper Functions ---

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
    """Generates the standard query filename based on parameters including data file."""
    return f"query-x{block_size}-n{number}-w{window}-v{value_range:.1f}-{boolean_logic.lower()}-s{selectivity}-{agg_type.upper()}.json"


# --- Core Functions ---

def parse_line(line: str) -> Optional[Dict[str, Any]]:
    """ Parses a single line from the data file. Returns a dict or None if parsing fails."""
    match = LINE_REGEX.match(line)
    if match:
        try:
            block_number = int(match.group(1))
            # Use float directly, Decimal was causing issues later with MongoDB types if not careful
            amount_str = match.group(2)
            amount = float(amount_str)
            addresses_str = match.group(3)
            # Filter out empty strings resulting from consecutive commas or trailing comma
            addresses = [addr.strip().strip("'\"") for addr in addresses_str.split(',') if addr.strip()]
            # Ensure no empty strings remain after stripping quotes
            addresses = [addr for addr in addresses if addr]
            return {"block_number": block_number, "amount": amount, "addresses": addresses}
        except ValueError as e:
            print(f"Warning: Data conversion error while parsing line '{line.strip()}': {e}")
            return None
        except Exception as e: # Catch broader exceptions during parsing
            print(f"Warning: Unexpected error while parsing line '{line.strip()}': {e}")
            return None
    else:
        # Only print warning if the line is not empty/whitespace
        if line.strip():
            print(f"Warning: Line format mismatch: '{line.strip()}'")
        return None

def load_data_to_mongo(file_path: str, db_uri: str, db_name: str) -> bool:
    """ Loads data from a file into a MongoDB collection (name derived from file_path), overwriting existing data. """
    client: Optional[pymongo.MongoClient] = None
    loaded_count = 0
    collection_name = get_collection_name(file_path) # Derive collection name

    try:
        if not os.path.exists(file_path):
            print(f"Error: Data file '{file_path}' not found. Please ensure the file exists.")
            return False

        client = pymongo.MongoClient(db_uri)
        db = client[db_name]
        collection = db[collection_name] # Use derived name

        # Clear existing data
        delete_result = collection.delete_many({})
        print(f"Cleared old data from collection '{db_name}.{collection_name}' ({delete_result.deleted_count} documents deleted).") # Use derived name

        documents: List[Dict[str, Any]] = []
        print(f"Starting to load data from '{file_path}' into collection '{collection_name}'...") # Use derived name
        with open(file_path, 'r') as f:
            for i, line in enumerate(f):
                parsed_doc = parse_line(line) # Use the refined parse_line
                if parsed_doc:
                    documents.append(parsed_doc)
                # Optional: Add logging for skipped lines here if needed, parse_line already prints warnings
                # else:
                #     print(f"  跳过第 {i+1} 行。") # parse_line handles warnings now

        if documents:
            print(f"Inserting {len(documents)} valid documents into '{collection_name}'...") # Use derived name
            insert_result = collection.insert_many(documents)
            loaded_count = len(insert_result.inserted_ids)
            print(f"Successfully loaded {loaded_count} documents into '{db_name}.{collection_name}'.") # Use derived name

            # Create indexes (consider checking if they exist first for efficiency)
            # It's often okay to run create_index multiple times; MongoDB handles it.
            print(f"Attempting to create/confirm indexes for collection '{collection_name}'...") # Use derived name
            try:
                # Add background=True for potentially large collections? Optional.
                collection.create_index([("block_number", pymongo.ASCENDING)], name="block_number_idx")
                collection.create_index([("amount", pymongo.ASCENDING)], name="amount_idx")
                collection.create_index([("addresses", pymongo.ASCENDING)], name="addresses_idx")
                print(f"Created/confirmed indexes for block_number, amount, addresses in '{collection_name}'.") # Use derived name
            except pymongo.errors.OperationFailure as e:
                # This might happen if indexes exist with different options, or due to permissions
                print(f"Error creating indexes (might already exist or permission issue): {e}")
        else:
            print(f"No valid data found in file '{file_path}' to load.")

        return True # Indicate success

    except pymongo.errors.ConnectionFailure as e:
        print(f"Database connection error: {e}")
        return False
    except FileNotFoundError: # Should be caught by os.path.exists, but good practice
        print(f"Error: Data file '{file_path}' not found.")
        return False
    except Exception as e:
        print(f"An unknown error occurred while loading data: {e}")
        return False
    finally:
        if client:
            client.close()
            # Only print connection closed if we actually loaded something or attempted to.
            if loaded_count > 0 or 'collection' in locals():
                 print(f"Database connection closed (operated on collection: {collection_name}).") # Use derived name


def query_mongo_data(db_uri: str, db_name: str, collection_name: str,
                     agg_type: str,
                     begin_block: Optional[int] = None, end_block: Optional[int] = None,
                     low_amount: Optional[float] = None, up_amount: Optional[float] = None,
                     address_filter: Optional[List[str]] = None,
                     logic_operator: str = "OR"
                     ) -> Optional[Union[int, float]]:
    """ Performs aggregation queries on the MongoDB collection. """
    result_value: Optional[Union[int, float]] = None
    client: Optional[pymongo.MongoClient] = None
    valid_agg_types = ["count", "sum", "max", "countdist"]
    agg_type_lower = agg_type.lower()

    if agg_type_lower not in valid_agg_types:
        print(f"Error: Invalid aggregation type: {agg_type}. Please use one of {valid_agg_types}.")
        return None # Return None for invalid agg_type

    try:
        client = pymongo.MongoClient(db_uri)
        db = client[db_name]
        collection = db[collection_name]

        # --- Build Match Stage ---
        match_conditions: Dict[str, Any] = {}

        # Block filter
        block_filter: Dict[str, int] = {}
        if begin_block is not None:
            block_filter["$gte"] = begin_block
        if end_block is not None:
            block_filter["$lte"] = end_block
        if block_filter:
            match_conditions["block_number"] = block_filter

        # Amount filter
        amount_filter: Dict[str, float] = {}
        # Ensure amounts are floats for MongoDB comparison if provided
        if low_amount is not None:
             try: amount_filter["$gte"] = float(low_amount)
             except (ValueError, TypeError): print(f"Warning: Invalid low_amount value {low_amount}, it will be ignored."); pass # Or return None/raise error
        if up_amount is not None:
             try: amount_filter["$lte"] = float(up_amount)
             except (ValueError, TypeError): print(f"Warning: Invalid up_amount value {up_amount}, it will be ignored."); pass # Or return None/raise error
        if amount_filter:
            match_conditions["amount"] = amount_filter

        # Address filter
        if address_filter and isinstance(address_filter, list) and len(address_filter) > 0:
            # Ensure all elements are strings
            clean_address_filter = [str(addr) for addr in address_filter if isinstance(addr, (str, int, float))] # Basic check
            op = logic_operator.upper()
            if op == "AND":
                # Remove duplicates for $all query
                unique_addresses = list(set(clean_address_filter))
                match_conditions["addresses"] = {"$all": unique_addresses}
            elif op == "OR":
                match_conditions["addresses"] = {"$in": clean_address_filter}
            else:
                print(f"Warning: Invalid address logical operator '{logic_operator}'. Address filter will be ignored.")
        # --- End Build Match Stage ---


        # --- Build Group Stage ---
        group_stage: Dict[str, Any] = {"_id": None} # Group all matched documents together
        if agg_type_lower == "count":
            group_stage["result"] = {"$sum": 1}
        elif agg_type_lower == "sum":
            # Ensure sum doesn't fail if amount is missing in some docs (though our parser should prevent this)
            group_stage["result"] = {"$sum": "$amount"}
        elif agg_type_lower == "max":
            group_stage["result"] = {"$max": "$amount"}
        elif agg_type_lower == "countdist":
            # Use $addToSet to collect unique amounts
            group_stage["distinct_amounts"] = {"$addToSet": "$amount"}
        # --- End Build Group Stage ---


        # --- Build Pipeline ---
        pipeline: List[Dict[str, Any]] = []
        if match_conditions: # Only add $match stage if there are conditions
            pipeline.append({"$match": match_conditions})
        pipeline.append({"$group": group_stage})

        # If countdist, add a stage to count the size of the set
        if agg_type_lower == "countdist":
            pipeline.append({
                "$project": {
                    "_id": 0,
                    # Calculate the size of the distinct_amounts array
                    # Use $ifNull to return 0 if distinct_amounts doesn't exist (no matching docs)
                    "result": {"$ifNull": [{"$size": "$distinct_amounts"}, 0]}
                }
            })
        # --- End Build Pipeline ---

        # --- Execute Aggregation ---
        aggregation_result = list(collection.aggregate(pipeline))

        # --- Process Result ---
        if aggregation_result:
            # Result is potentially an int (count) or float (sum, max, min)
            result_value = aggregation_result[0].get("result")
        else:
            # If no documents match, return appropriate default
            if agg_type_lower in ["count", "sum"]:
                result_value = 0
            else: # For min/max, None is more appropriate if no docs match
                result_value = None

    except pymongo.errors.ConnectionFailure as e:
        print(f"Database connection error: {e}")
        # Optionally, re-raise or return a specific error code/value
    except ValueError as e: # Catch potential value errors during type conversions etc.
        print(f"Value error during query processing: {e}")
    except Exception as e:
        print(f"An unknown error occurred while querying data: {e}")
    finally:
        if client:
            client.close()
            # print("Database connection closed after query.") # Optional: less verbose

    return result_value

def delete_database(db_uri: str, db_name: str) -> None:
    """ Deletes the specified MongoDB database. """
    client: Optional[pymongo.MongoClient] = None
    try:
        client = pymongo.MongoClient(db_uri)
        db_list = client.list_database_names()
        if db_name in db_list:
            client.drop_database(db_name)
            print(f"Database '{db_name}' deleted successfully.")
        else:
            print(f"Database '{db_name}' does not exist, no need to delete.")
    except pymongo.errors.ConnectionFailure as e:
        print(f"Database connection error: {e}")
    except Exception as e:
        print(f"Error deleting database '{db_name}': {e}")
    finally:
        if client:
            client.close()
            print("Database connection closed.")


# --- Query File Execution Helper Functions ---

def _parse_range(query_spec: Dict[str, Any], query_num: int) -> Tuple[Optional[float], Optional[float]]:
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
                except (ValueError, TypeError):
                    print(f"Warning: Query {query_num} - Invalid range value (must be number or None): {range_spec[0]}. Skipping range filter.")
                    low_amount, up_amount = None, None # Reset on error
            else:
                 print(f"Warning: Query {query_num} - Invalid 'range' inner list format (needs [low, high]): {range_spec[0]}. Skipping range filter.")
        else:
            print(f"Warning: Query {query_num} - Invalid 'range' format (needs list): {range_spec}. Skipping range filter.")
    return low_amount, up_amount


def _parse_keyword_exp(query_spec: Dict[str, Any], query_num: int) -> Tuple[str, List[str]]:
    """Parses the 'keyword_exp' field, returning logic operator and address list."""
    logic_operator = "OR"  # Default
    address_list_for_query: List[str] = []
    keyword_exp = query_spec.get("keyword_exp")

    if keyword_exp and isinstance(keyword_exp, dict):
        address_inputs: Optional[List[Any]] = None
        found_logic = False

        # Check for 'and' or 'or' keys, case-insensitive check might be safer if needed
        if "and" in keyword_exp and isinstance(keyword_exp["and"], list):
            logic_operator, address_inputs, found_logic = "AND", keyword_exp["and"], True
        elif "or" in keyword_exp and isinstance(keyword_exp["or"], list):
            logic_operator, address_inputs, found_logic = "OR", keyword_exp["or"], True
        else:
            print(f"Warning: Query {query_num} - 'keyword_exp' must contain 'and' or 'or' key with a list value. Skipping keyword filter.")

        if found_logic and address_inputs is not None:
            extracted_addresses: List[str] = []
            for item in address_inputs:
                # Expecting {"input": "'address_value'"}
                if isinstance(item, dict) and "input" in item and isinstance(item["input"], str):
                    addr = item["input"].strip().strip("'\"") # Strip quotes and whitespace
                    if addr: # Ensure address is not empty after stripping
                        extracted_addresses.append(addr)
                    else:
                        print(f"Warning: Query {query_num} - Empty or invalid 'input' value in 'keyword_exp': {item}")
                else:
                    print(f"Warning: Query {query_num} - Invalid item format in 'keyword_exp' (needs {{'input': 'addr'}}): {item}")

            if extracted_addresses:
                address_list_for_query = extracted_addresses
            elif found_logic: # Only warn if we expected addresses but got none
                 print(f"Warning: Query {query_num} - Failed to extract valid addresses from 'keyword_exp'.")

    return logic_operator, address_list_for_query


def execute_queries_from_file(query_file: str, db_uri: str, db_name: str) -> float:
    """ Reads query definitions from a JSON file and executes them against the collection derived from the filename. """
    try:
        if not os.path.exists(query_file):
            print(f"Error: Query file '{query_file}' not found.")
            return
        with open(query_file, 'r') as f:
            queries_data = json.load(f)
    except json.JSONDecodeError as e:
        print(f"Error: Failed to parse query file '{query_file}': {e}")
        return
    except IOError as e:
        print(f"Error: Error reading query file '{query_file}': {e}")
        return
    except Exception as e: # Catch other potential errors during file read/load
        print(f"An unknown error occurred while reading or parsing the query file: {e}")
        return

    # --- Derive collection name from query filename ---
    base_query_filename = os.path.basename(query_file)
    # Expecting format: query-{data_id}-n{number}-...json
    match = re.match(r"query-(.*?)-n\d+.*\.json", base_query_filename)
    if not match:
        print(f"Error: Cannot parse data identifier from query filename '{base_query_filename}'.")
        print(f"       Expected format: 'query-DATAID-nNUMBER-...' (e.g., 'query-eth_32-n100-...')")
        return
    data_id = match.group(1)
    # Reconstruct a dummy path for the helper to get the collection name
    collection_name = get_collection_name(f"{data_id}.dat") # Assume .dat extension was original
    print(f"\n--- Executing queries from '{query_file}' (Target collection: '{collection_name}') ---")
    # --- End Derive collection name ---


    if not isinstance(queries_data, list):
        print(f"Error: The top-level structure of query file '{query_file}' must be a JSON array/list.")
        return

    query_execution_started = False
    total_duration_ms = 0.0 # Initialize total duration
    executed_query_count = 0 # Initialize count of executed queries

    for i, query_spec in enumerate(queries_data):
        query_num = i + 1
        print(f"\n--- Query {query_num} ---")
        if not isinstance(query_spec, dict):
            print(f"Warning: Skipping invalid query entry {query_num} (not a dictionary): {query_spec}")
            continue

        query_execution_started = True # Mark that we started processing at least one query spec

        try:
            # --- Extract common parameters ---
            agg_type = query_spec.get("agg_type")
            if not agg_type or not isinstance(agg_type, str):
                print(f"Warning: Skipping query {query_num} - Missing or invalid 'agg_type' field.")
                continue

            begin_block = query_spec.get("start_blk")
            end_block = query_spec.get("end_blk")

            # --- Parse range and keyword_exp using helper functions ---
            low_amount, up_amount = _parse_range(query_spec, query_num)
            logic_operator, address_list_for_query = _parse_keyword_exp(query_spec, query_num)

            # --- Print query details ---
            print(f"  Aggregation Type: {agg_type.upper()}")
            print(f"  Block Range: {begin_block if begin_block is not None else 'N/A'} - {end_block if end_block is not None else 'N/A'}")
            print(f"  Amount Range: {low_amount if low_amount is not None else 'N/A'} - {up_amount if up_amount is not None else 'N/A'}")
            print(f"  Address Logic: {logic_operator}")
            print(f"  Address Filter: {address_list_for_query if address_list_for_query else 'None'}")

            # --- Execute the query ---
            start_time = time.perf_counter() # Record start time
            result = query_mongo_data(
                db_uri, db_name, collection_name, # Pass derived collection_name
                agg_type=agg_type,
                begin_block=begin_block,
                end_block=end_block,
                low_amount=low_amount,
                up_amount=up_amount,
                address_filter=address_list_for_query,
                logic_operator=logic_operator
            )
            end_time = time.perf_counter() # Record end time
            duration_ms = (end_time - start_time) * 1000 # Calculate duration in milliseconds

            # --- Add to total duration and increment count ---
            total_duration_ms += duration_ms
            executed_query_count += 1

            # --- Print result ---
            # Format float results nicely
            if isinstance(result, float):
                print(f"Result: {result:.4f}")
            else:
                print(f"Result: {result if result is not None else 'No matching documents'}") # Handle None result explicitly
            print(f"Execution Time: {duration_ms:.2f} ms") # Print execution time

        except Exception as e:
            # Catch unexpected errors during the processing of a single query
            print(f"Unexpected error processing query {query_num}: {query_spec}\nError details: {e}")
            # Continue to the next query

    if query_execution_started:
        print("\n--- JSON Query Execution Finished ---")
        # --- Calculate and print average time ---
        if executed_query_count > 0:
            average_time_ms = total_duration_ms / executed_query_count
            print(f"Executed {executed_query_count} queries, Average Execution Time: {average_time_ms:.2f} ms")
        else:
            print("No queries were executed successfully, cannot calculate average time.") # Handle case where no queries were executed successfully
    elif isinstance(queries_data, list): # Check if it was a list but maybe empty
         print("Query file is empty or all entries were invalid.")
    
    return average_time_ms


def get_data_bounds(db_uri: str, db_name: str, collection_name: str) -> Optional[Tuple[int, int, float, float]]:
    """ Retrieves min/max block number and min/max amount from the specified collection. """
    client: Optional[pymongo.MongoClient] = None
    try:
        client = pymongo.MongoClient(db_uri)
        db = client[db_name]
        collection = db[collection_name] # Use specified collection name

        # Check if collection is empty first, as aggregate might return weird results or error
        count = collection.count_documents({})
        if count == 0:
            print(f"Info: Collection '{collection_name}' is empty, cannot get data bounds. Please load data first.")
            return None # Return None, it's not strictly an error, just no data

        # Pipeline to get all bounds in one go
        pipeline = [
            {'$group': {
                '_id': None, # Group all documents
                'min_blk': {'$min': '$block_number'},
                'max_blk': {'$max': '$block_number'},
                'min_amt': {'$min': '$amount'},
                'max_amt': {'$max': '$amount'}
            }}
        ]
        bounds_result = list(collection.aggregate(pipeline))

        # Check if aggregation returned any result (it should if collection is not empty)
        if not bounds_result: # If the list is empty, aggregation failed unexpectedly
            print("Error: Could not get data bounds via aggregation (query returned no result even though collection is not empty).")
            return None

        # Now we know bounds_result[0] exists
        bounds = bounds_result[0]
        min_blk = bounds.get('min_blk')
        max_blk = bounds.get('max_blk')
        min_amt = bounds.get('min_amt')
        max_amt = bounds.get('max_amt')

        # Validate that we got numbers back
        if not all(isinstance(v, (int, float)) for v in [min_blk, max_blk, min_amt, max_amt]):
             print(f"Error: Data bounds query returned invalid or missing values: Blk({min_blk}, {max_blk}), Amt({min_amt}, {max_amt})")
             return None

        # Ensure correct types (MongoDB might return int for amounts if all are whole numbers)
        min_blk, max_blk = int(min_blk), int(max_blk)
        min_amt, max_amt = float(min_amt), float(max_amt)


        print(f"Data Bounds: Block [{min_blk}, {max_blk}], Amount [{min_amt:.4f}, {max_amt:.4f}]")
        return min_blk, max_blk, min_amt, max_amt

    except pymongo.errors.ConnectionFailure as e:
        print(f"Database connection error: {e}")
        return None
    except Exception as e:
        print(f"An unknown error occurred while getting data bounds for collection '{collection_name}': {e}") # Include collection name in error
        return None
    finally:
        if client:
            client.close()

# --- Query Generation Helper Functions ---



# --- Helper function for nested OR keyword expressions ---
def _build_nested_or_exp(addresses: List[str]) -> Dict[str, Any]:
    """Recursively builds a nested OR expression for the query file."""
    if not addresses:
        return {} # Return empty dict if no addresses
    
    # Format address helper
    def format_input(addr: str) -> Dict[str, str]:
        return {"input": f"'{addr}'"}

    count = len(addresses)
    if count == 1:
        return format_input(addresses[0])
    elif count == 2:
        return {"or": [format_input(addresses[0]), format_input(addresses[1])]}
    else:
        # Recursive step: take the first, nest the rest
        return {"or": [format_input(addresses[0]), _build_nested_or_exp(addresses[1:])]}

# Updated: generate_query_file uses data_file_path and derived collection name
def generate_query_file(filename: str, # Output query filename
                        data_file_path: str, # Input data file to base queries on
                        number: int, window: int, value_range_ratio: float,
                        boolean_logic: str, selectivity: int, agg_type: str,
                        db_uri: str, db_name: str) -> None:
    """
    Generates a JSON query file with n random queries using a single specified aggregation type
    and window-specific keyword selection, based on the data in the collection derived from data_file_path.
    """
    collection_name = get_collection_name(data_file_path) # Derive collection name
    print(f"\n--- Starting query file generation: {filename} (Based on data file: {data_file_path}, Target collection: {collection_name}) ---")
    print(f"Parameters: N={number}, W={window}, V={value_range_ratio:.2f}, Logic={boolean_logic.upper()}, S={selectivity}, Agg={agg_type.upper()}")

    # 1. Get Data Bounds from the correct collection
    bounds_data = get_data_bounds(db_uri, db_name, collection_name) # Use derived name
    if bounds_data is None:
        print(f"Error: Cannot get data bounds from collection '{collection_name}', generation aborted. Ensure data is loaded into this collection.")
        return
    min_blk, max_blk, min_amt, max_amt = bounds_data
    print(f"Data Bounds: Block [{min_blk}, {max_blk}], Amount [{min_amt:.4f}, {max_amt:.4f}]")

    target_amt_width = max_amt * value_range_ratio

    # 3. Connect to DB and get the specific collection
    client: Optional[pymongo.MongoClient] = None
    collection: Optional[Collection] = None
    try:
        client = pymongo.MongoClient(db_uri)
        db = client[db_name]
        collection = db[collection_name] # Get the specific collection object
        # Test connection early
        client.admin.command('ping')
    except pymongo.errors.ConnectionFailure as e:
        print(f"Database connection error: {e}")
        if client: client.close()
        return
    except Exception as e: # Catch other potential DB connection errors
        print(f"Error connecting to database or accessing collection '{collection_name}': {e}") # Include collection name
        if client: client.close()
        return


    # 4. Generate Queries
    queries_list: List[Dict[str, Any]] = []
    print(f"Generating {number} queries...")

    # Ensure collection is not None before proceeding (should be caught by try/except above, but defensive check)
    if collection is None:
        print(f"Error: Failed to get collection '{collection_name}' object.")
        if client: client.close()
        return

    for i in range(number):
        # --- Determine Block Range ---
        # Ensure start_blk allows for the full window width within the available range
        max_start_blk = max(min_blk, max_blk - window) # Ensure start_blk is at least min_blk
        start_blk = random.randint(min_blk, max_start_blk)
        assert(max_start_blk >= min_blk)
        end_blk = start_blk + window - 1
        assert(start_blk < end_blk)

        # --- Determine Amount Range ---
        low_amount, up_amount = int(min_amt), int(min_amt) # Default if range is zero
        max_start_amount = max(min_amt, max_amt - target_amt_width) # Ensure low_amount is at least min_amt
        # Convert random float to integer
        low_amount = int(random.uniform(min_amt, max_start_amount))
        up_amount = int(low_amount + target_amt_width)
        assert(up_amount <= max_amt)

                # --- Select Keywords (Addresses) using new logic ---
        selected_addresses: List[str] = []
        k = selectivity # Target number of unique addresses

        while len(selected_addresses) < k:
            rand_blk = random.randint(min_blk, max_blk)

            # 2. Find transactions in that block (fetch only addresses)
            # Limit the number fetched for efficiency if blocks can be huge
            cursor = collection.find(
                {"block_number": rand_blk},
                {"addresses": 1, "_id": 0}
            )
            block_transactions = list(cursor) # Evaluate cursor to list

            # 3. Select random transaction
            chosen_tx = random.choice(block_transactions)
            tx_addresses = chosen_tx.get("addresses")

            # 4. Select random address from the transaction
            addr = random.choice(tx_addresses)
            selected_addresses.append(addr)
        # --- End Address Selection ---


        # --- Format keyword_exp based on selected_addresses ---
        keyword_exp: Dict[str, Any] = {}
        if selected_addresses: # Proceed only if we found at least one address
            logic_key = boolean_logic.lower()
            current_k = len(selected_addresses) # The actual number of unique addresses found

            # Define helper inside or ensure it's available
            def format_input(addr: str) -> Dict[str, str]:
                 return {"input": f"'{addr}'"}

            if logic_key == "or":
                if current_k > 1:
                    # Use the recursive helper for nested OR structure
                    keyword_exp = _build_nested_or_exp(selected_addresses)
                elif current_k == 1:
                    # Directly format the single address
                    keyword_exp = format_input(selected_addresses[0])
                # If current_k is 0, keyword_exp remains {}

            elif logic_key == "and":
                 # AND logic uses a flat list
                 # If current_k >= 1: (handles case where k=0 or no addresses found)
                 if current_k > 0:
                    keyword_exp = {
                        "and": [format_input(addr) for addr in selected_addresses]
                    }
                 # If current_k is 0, keyword_exp remains {}

            # If logic_key is somehow neither 'or' nor 'and', keyword_exp remains {}


        # --- Construct Query Object ---
        query_obj = {
            "agg_type": agg_type.upper(), # Standardize to uppercase
            "start_blk": start_blk,
            "end_blk": end_blk,
            # Use integer values, no rounding needed
            "range": [[low_amount, up_amount]],
        }

        # Conditionally add keyword_exp if selectivity is non-zero
        if selectivity != 0:
            query_obj["keyword_exp"] = keyword_exp

        queries_list.append(query_obj)

    # 5. Close DB Connection
    if client:
        client.close()

    # 6. Write to JSON File
    try:
        with open(filename, 'w') as f:
            json.dump(queries_list, f, indent=2)
        print(f"Successfully wrote {number} queries to file: '{filename}'")
    except IOError as e:
        print(f"Error: Failed to write query file '{filename}': {e}")
    except Exception as e:
        print(f"An unknown error occurred while writing query file: {e}")

# --- Main execution block ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Load, query blockchain transaction data, generate query files, or delete the database using MongoDB. Collection names and data file names are dynamically generated based on block size.",
        formatter_class=argparse.RawTextHelpFormatter # Keep help formatting
    )

    subparsers = parser.add_subparsers(dest='action', required=True, help='Select the action to perform (load, query, create-query, delete)')

    # --- Load Action ---
    parser_load = subparsers.add_parser('load', help="Load data from eth-{block_size}.dat based on block size (will clear or create the corresponding collection)")
    parser_load.add_argument('-x', '--block-size', type=lambda v: _validate_positive_int(v, '-x/--block-size'), required=True,
                             help="Specify the block size of the data file (e.g., 32 loads eth-32.dat). Collection name will be generated based on this.")

    # --- Query Action ---
    parser_query = subparsers.add_parser('query',
                                         help=f"""Execute queries from the specified JSON file.
                                        The query file itself should contain a data source identifier (e.g., query-eth_32-n100...)
                                        The script will infer the MongoDB collection to query based on this.""",
                                       formatter_class=argparse.RawTextHelpFormatter)
    parser_query.add_argument('-f', '--file', type=str, required=False,
                             help='Specify the query JSON file to execute. The filename needs to follow the `query-DATAID-nNUMBER...` format.')
    parser_query.add_argument('-c', '--case', type=lambda v: _validate_positive_int(v, '-c/--case'), default=1,
                             help="Case 0: Single query using provided parameters.\n"
                                  "Case 1: Vary time window, testing both OR & AND logic.\n"
                                  "Case 2: Vary value range.\n"
                                  "Case 3: Vary keyword number (selectivity).\n"
                                  "Case 4: Vary block size (data source).\n"
                                  "Case 5: Vary aggregation type.")

    # --- Delete Action ---
    parser_delete = subparsers.add_parser('delete', help=f"Delete the entire database '{DB_NAME}' (including all collections)")
    # No specific arguments for delete

    # --- Create Query File Action ---
    parser_create = subparsers.add_parser('create-query', help="Generate a query file for data with a specific block size (eth-{block_size}.dat)")
    parser_create.add_argument('-x', '--block-size', type=lambda v: _validate_positive_int(v, '-x/--block-size'), default=128,
                               help="Specify the block size of the base data for query generation (e.g., 128)." )
    parser_create.add_argument('-n', '--number', type=lambda v: _validate_positive_int(v, '-n/--number'), default=100,
                               help="Number of queries to generate (default: 100)")
    parser_create.add_argument('-w', '--window', type=lambda v: _validate_non_negative_int(v, '-w/--window'), default=800,
                               help="Time window size (default: 800)")
    parser_create.add_argument('-v', '--value-range', type=lambda v: _validate_range_float(v, '-v/--value-range'), default=0.5,
                               help="Amount range ratio [0.0-1.0] (default: 0.5)")
    parser_create.add_argument('-b', '--boolean', choices=['AND', 'OR'], default='OR', type=str.upper, # Ensure uppercase
                               help="Boolean logic for keyword queries (default: OR)")
    parser_create.add_argument('-s', '--selectivity', type=lambda v: _validate_positive_int(v, '-s/--selectivity'), default=2,
                               help="Number of keywords to select (selectivity) (default: 2)")
    parser_create.add_argument('-a', '--agg-type', choices=['MAX', 'COUNT', 'SUM', 'COUNTDIST'], default='COUNT', type=str.upper, # Ensure uppercase
                               help="Aggregation type (default: COUNT)")
    parser_create.add_argument('-c', '--case', type=lambda v: _validate_positive_int(v, '-c/--case'), default=1,
                               help="Case 0: Generate a single query file using provided parameters.\n"
                                    "Case 1: Generate files varying time window, for both OR & AND logic.\n"
                                    "Case 2: Generate files varying value range.\n"
                                    "Case 3: Generate files varying keyword number (selectivity).\n"
                                    "Case 4: Generate files varying block size (data source).\n"
                                    "Case 5: Generate files varying aggregation type.")

    # --- Parse Arguments ---
    try:
        args = parser.parse_args()
    except argparse.ArgumentTypeError as e:
         # Exit gracefully if validation helpers raise ArgumentTypeError
         parser.error(str(e))


    # --- Execute Action ---
    if args.action == 'load':
        # Construct data file path from block size
        data_file_path = f"../data/eth-{args.block_size}.dat"
        print(f"\n--- Selected Action: Load Data (from {data_file_path}) ---")
        # Pass the constructed path
        success = load_data_to_mongo(data_file_path, MONGO_URI, DB_NAME)
        print(f"--- Data Loading Action {'Completed' if success else 'Failed'} ---")

    elif args.action == 'query':
        print("\n--- Selected Action: Execute Queries ---")
        
        block_size = 128
        number = 100
        window = 800
        value_range = 0.5
        boolean = 'OR'
        selectivity = 2
        agg_type = 'COUNT'

        if args.case == 0:
            query_file_to_use = args.file # Already required by argparse
            print(f"Using query file: '{query_file_to_use}'")
            # Basic check if file exists before proceeding
            if not os.path.exists(query_file_to_use):
                print(f"Error: Specified query file '{query_file_to_use}' does not exist.")
            else:
                # Pass only query file, db_uri, db_name. Collection name derived inside.
                execute_queries_from_file(query_file_to_use, MONGO_URI, DB_NAME)
        elif args.case == 1:
            booleans = ['AND', 'OR']
            windows = [400, 800, 1200, 1600, 2000]

            for b in booleans:
                output_file = f'res/window_{b.lower()}.csv'
                output_dir = os.path.dirname(output_file)
                os.makedirs(output_dir, exist_ok=True)
                # Write header to output file
                with open(output_file, 'w') as f:
                    f.write("time window,\tquery time\\n")

                    for w in windows:
                        query_file = f'query/window_{b.lower()}/' + _generate_query_filename(
                            block_size, number, w, value_range, b, selectivity, agg_type
                        )
                        query_time = execute_queries_from_file(query_file, MONGO_URI, DB_NAME)
                        f.write(f"{w},\t{query_time}\\n")
        elif args.case == 2:
            value_ranges = [0.1, 0.3, 0.5, 0.7, 0.9]
            selectivity = 0

            output_file = f'res/value_range.csv'
            output_dir = os.path.dirname(output_file)
            os.makedirs(output_dir, exist_ok=True)
            # Write header to output file
            with open(output_file, 'w') as f:
                f.write("value range,\tquery time\\n")

                for v in value_ranges:
                    query_file = f'query/value_range/' + _generate_query_filename(
                        block_size, number, window, v, boolean, selectivity, agg_type
                    )
                    query_time = execute_queries_from_file(query_file, MONGO_URI, DB_NAME)
                    f.write(f"{v},\t{query_time}\\n")
        elif args.case == 3:
            keywords = [2, 4, 8, 16, 32]

            output_file = f'res/keyword.csv'
            output_dir = os.path.dirname(output_file)
            os.makedirs(output_dir, exist_ok=True)
            # Write header to output file
            with open(output_file, 'w') as f:
                f.write("keyword number,\tquery time\\n")

                for k in keywords:
                    query_file = f'query/keyword/' + _generate_query_filename(
                        block_size, number, window, value_range, boolean, k, agg_type
                    )
                    query_time = execute_queries_from_file(query_file, MONGO_URI, DB_NAME)
                    f.write(f"{k},\t{query_time}\\n")
        elif args.case == 4:
            block_sizes = [32, 64, 128, 256, 512]

            output_file = f'res/blk_size.csv'
            output_dir = os.path.dirname(output_file)
            os.makedirs(output_dir, exist_ok=True)
            #  Write header to output file
            with open(output_file, 'w') as f:
                f.write("blk size,\tquery time\\n")

                for b in block_sizes:
                    query_file = f'query/blk_size/' + _generate_query_filename(
                        b, number, window, value_range, boolean, selectivity, agg_type
                    )
                    query_time = execute_queries_from_file(query_file, MONGO_URI, DB_NAME)
                    f.write(f"{b},\t{query_time}\\n")
        elif args.case == 5:
            agg_types = ['MAX', 'COUNT', 'SUM', 'COUNTDIST']

            output_file = f'res/agg_type.csv'
            output_dir = os.path.dirname(output_file)
            os.makedirs(output_dir, exist_ok=True)
            #  Write header to output file
            with open(output_file, 'w') as f:
                f.write("agg type,\tquery time\\n")

                for a in agg_types:
                    query_file = f'query/agg_type/' + _generate_query_filename(
                        block_size, number, window, value_range, boolean, selectivity, a
                    )
                    query_time = execute_queries_from_file(query_file, MONGO_URI, DB_NAME)
                    f.write(f"{a},\t{query_time}\\n")

    elif args.action == 'delete':
        print("\n--- Selected Action: Delete Database ---")
        # Add confirmation prompt? (Optional but recommended for destructive actions)
        confirm = input(f"Are you sure you want to delete the entire database '{DB_NAME}' (including all dynamically generated collections)? (yes/no): ").lower()
        if confirm == 'yes':
            delete_database(MONGO_URI, DB_NAME)
            print(f"--- Database deletion attempt finished ---")
        else:
            print("Database deletion cancelled.")

    elif args.action == 'create-query':
        print("\n--- Selected Action: Create Query File ---")

        base_dir = "query"

        # Replace switch with if/elif
        if args.case == 0:
            # Generate the output filename using the *constructed* data file path
            output_filename = f'{base_dir}/{args.block_size}/' + _generate_query_filename(
                args.block_size, args.number, args.window, args.value_range, args.boolean, args.selectivity, args.agg_type
            )
            # Ensure the output directory exists
            output_dir = os.path.dirname(output_filename)
            os.makedirs(output_dir, exist_ok=True)
            print(f"Generated query filename will be: '{output_filename}' (Based on data: eth-{args.block_size}.dat)")
            generate_query_file(
                filename=output_filename, # The generated output filename
                data_file_path= f"eth-{args.block_size}.dat",
                number=args.number,
                window=args.window,
                value_range_ratio=args.value_range,
                boolean_logic=args.boolean,
                selectivity=args.selectivity,
                agg_type=args.agg_type.upper(), # Pass uppercase agg_type
                db_uri=MONGO_URI,
                db_name=DB_NAME
            )

        elif args.case == 1:
            booleans = ['and', 'or']
            window = [400, 800, 1200, 1600, 2000]
            for boolean in booleans:
                for w in window:
                    output_filename = f'{base_dir}/window_{boolean}/' + _generate_query_filename(
                        args.block_size, args.number, w, args.value_range, boolean, args.selectivity, args.agg_type
                    )
                    # Ensure the output directory exists
                    output_dir = os.path.dirname(output_filename)
                    os.makedirs(output_dir, exist_ok=True)
                    print(f"Generated query filename will be: '{output_filename}' (Based on data: eth-{args.block_size}.dat)")
                    generate_query_file(
                        filename=output_filename, # The generated output filename
                        data_file_path= f"eth-{args.block_size}.dat",
                        number=args.number,
                        window=w,
                        value_range_ratio=args.value_range,
                        boolean_logic=boolean,
                        selectivity=args.selectivity,
                        agg_type=args.agg_type.upper(), # Pass uppercase agg_type
                        db_uri=MONGO_URI,
                        db_name=DB_NAME
                    )
        
        elif args.case == 2:
            value_ranges = [0.1, 0.3, 0.5, 0.7, 0.9]
            for value_range in value_ranges:
                output_filename = f'{base_dir}/value_range/' + _generate_query_filename(
                    args.block_size, args.number, args.window, value_range, args.boolean, args.selectivity, args.agg_type
                )
                # Ensure the output directory exists
                output_dir = os.path.dirname(output_filename)
                os.makedirs(output_dir, exist_ok=True)
                print(f"Generated query filename will be: '{output_filename}' (Based on data: eth-{args.block_size}.dat)")
                generate_query_file(
                    filename=output_filename, # The generated output filename
                    data_file_path= f"eth-{args.block_size}.dat",
                    number=args.number,
                    window=args.window,
                    value_range_ratio=value_range,
                    boolean_logic=args.boolean,
                    selectivity=args.selectivity,
                    agg_type=args.agg_type.upper(), # Pass uppercase agg_type
                    db_uri=MONGO_URI,
                    db_name=DB_NAME
                )        

        elif args.case == 3:
            keywords = [1, 2, 4, 8, 16]
            for keyword in keywords:
                    output_filename = f'{base_dir}/keyword/' + _generate_query_filename(
                        args.block_size, args.number, args.window, args.value_range, args.boolean, keyword, args.agg_type
                    )
                    # Ensure the output directory exists
                    output_dir = os.path.dirname(output_filename)
                    os.makedirs(output_dir, exist_ok=True)
                    print(f"Generated query filename will be: '{output_filename}' (Based on data: eth-{args.block_size}.dat)")
                    generate_query_file(
                        filename=output_filename, # The generated output filename
                        data_file_path= f"eth-{args.block_size}.dat",
                        number=args.number,
                        window=args.window,
                        value_range_ratio=args.value_range,
                        boolean_logic=args.boolean,
                        selectivity=keyword,
                        agg_type=args.agg_type.upper(), # Pass uppercase agg_type
                        db_uri=MONGO_URI,
                        db_name=DB_NAME
                    )     

        elif args.case == 4:
            block_sizes = [32, 64, 128, 256, 512]
            for block_size in block_sizes:
                output_filename = f'{base_dir}/blk_size/' + _generate_query_filename(
                    block_size, args.number, args.window, args.value_range, args.boolean, args.selectivity, args.agg_type
                )
                # Ensure the output directory exists
                output_dir = os.path.dirname(output_filename)
                os.makedirs(output_dir, exist_ok=True)
                print(f"Generated query filename will be: '{output_filename}' (Based on data: eth-{args.block_size}.dat)")
                generate_query_file(
                    filename=output_filename, # The generated output filename
                    data_file_path= f"eth-{block_size}.dat",
                    number=args.number,
                    window=args.window,
                    value_range_ratio=args.value_range,
                    boolean_logic=args.boolean,
                    selectivity=args.selectivity,
                    agg_type=args.agg_type.upper(), # Pass uppercase agg_type
                    db_uri=MONGO_URI,
                    db_name=DB_NAME
                )
        
        elif args.case == 5:
            agg_types = ['MAX', 'COUNT', 'SUM', 'COUNTDIST']
            for agg_type in agg_types:
                output_filename = f'{base_dir}/agg_type/' + _generate_query_filename(
                    args.block_size, args.number, args.window, args.value_range, args.boolean, args.selectivity, agg_type
                )
                # Ensure the output directory exists
                output_dir = os.path.dirname(output_filename)
                os.makedirs(output_dir, exist_ok=True)
                print(f"Generated query filename will be: '{output_filename}' (Based on data: eth-{args.block_size}.dat)")
                generate_query_file(
                    filename=output_filename, # The generated output filename
                    data_file_path= f"eth-{args.block_size}.dat",
                    number=args.number,
                    window=args.window,
                    value_range_ratio=args.value_range,
                    boolean_logic=args.boolean,
                    selectivity=args.selectivity,
                    agg_type=agg_type.upper(), # Pass uppercase agg_type
                    db_uri=MONGO_URI,
                    db_name=DB_NAME
                )

        print("--- Create Query File Action Completed ---")

    else:
        # This should not be reachable due to 'required=True' in subparsers
        print("Error: Unknown action.")
        parser.print_help()