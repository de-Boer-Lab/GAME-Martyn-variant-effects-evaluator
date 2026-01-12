'''Handle Loading and Validating Evaluator Input/Request Data'''

import os
import json
import msgpack
from collections import Counter
import functools
import pandas as pd
from config import EVALUATOR_INPUT_PATH

class DuplicateKeysError(ValueError):
    """Raised when duplicate keys are found in a JSON object."""
    pass

# Internal helper function to detect duplicates during JSON parsing
def _detect_duplicates(pairs, duplicate_keys_state):

    """
    Detects duplicate keys during JSON parsing and counts occurrences of each key.

    This function intercepts the key-value pairs provided by `json.loads` and ensures that
    duplicate keys are flagged. It constructs the dictionary normally but counts how often
    each key appears, recording any keys that occur more than once.

    Args:
        pairs (list of tuple): A list of key-value pairs at the current level of the JSON.
        duplicate_keys_state (dict): The dictionary to update with any duplicates found.

    Returns:
        result_dict: A dictionary created from the key-value pairs.
    """

    # Use a local Counter to count occurrences of keys at this level
    local_counts = Counter()
    result_dict = {}
    for key, value in pairs:
        # Increment the count for each key
        local_counts[key] += 1
        # If the key is a duplicate, record it in the duplicate_keys dictionary
        if local_counts[key] > 1:
            duplicate_keys_state[key] = local_counts[key]
        # Add the key-value pair to the resulting dictionary
        result_dict[key] = value
    return result_dict

def _process_results(data, duplicate_keys):
    """
    Checks the duplicate_keys dictionary and prints a report.

    Args:
        data (dict): The dictionary of parsed data. 
        duplicate_keys (dict): The dictionary of duplicates.

    Returns:
        data or None: The parsed data if no duplicates. None, if duplicates are found.
    """
    # Report duplicates if any were found
    if duplicate_keys:
        print("Duplicate keys found:")
        error_messages = [f"Key: '{key}', Count: {count}" for key, count in duplicate_keys.items()]
        raise DuplicateKeysError(f"Duplicate keys found:\n" + "\n".join(error_messages))
    else:
        print("No duplicates found.")
        return data # Return the parsed data if no duplicates.


# Function to check for duplicate keys in JSON object

def check_duplicates_from_string(json_string):

    """
    Parses a JSON string to detect and report any duplicate keys at the same level in the same object.
    This function ensures that no keys are silently overwritten in dictionaries.

    The function uses a helper to track the number of times each key appears during parsing,
    leveraging the `object_pairs_hook` parameter of `json.loads()` to intercept key-value pairs
    before they are processed into a dictionary. If duplicates are detected at any level, they
    are reported with their counts. Keys reused in separate objects within arrays (e.g. lists) 
    are not considered duplicates.

    Args:
        json_string (str): The JSON content as a string to parse and check for duplicates.

    Raises:
        json.JSONDecodeError: If the string is not valid JSON.
        DuplicateKeysError: If duplicate keys are found in the JSON structure.

    Returns:
        dict: The parsed data if no errors or duplicates are found.
    """

    # Initialize a dictionary to track duplicate keys and their counts
    duplicate_keys = {}
    
    # Create a 1-argument hook callable by "freezing" the duplicate_keys dict
    # as the second argument to the helper.
    hook = functools.partial(_detect_duplicates, duplicate_keys_state=duplicate_keys)

    # Parse the JSON string using the helper to track duplicates
    data = json.loads(json_string, object_pairs_hook=hook)
    
    return _process_results(data, duplicate_keys)
    
# Function for check for duplicate keys if input file is in JSON format

def check_duplicates_from_json(json_file_path):
    """
    Parses a JSON file to detect and report any duplicate keys at the same level in the same object.
    This function ensures that no keys are silently overwritten in dictionaries.

    The function uses a helper to track the number of times each key appears during parsing,
    leveraging the `object_pairs_hook` parameter of `json.load()` to intercept key-value pairs 
    before they are processed into a dictionary. If duplicates are detected at any level, they
    are reported with their counts and paths. Keys reused in separate objects within arrays 
    (e.g. lists) are not considered duplicates.

    Raises:
        FileNotFoundError: If the specified file does not exist.
        json.JSONDecodeError: If the file content is not valid JSON.
        DuplicateKeysError: If duplicate keys are found in the JSON structure.

    Returns:
        dict: The parsed data if no errors or duplicates are found.
    """

    # Initialize a dictionary to track duplicate keys and their counts
    duplicate_keys = {}
    
    # Create a 1-argument hook callable by "freezing" the duplicate_keys dict
    # as the second argument to the helper.
    hook = functools.partial(_detect_duplicates, duplicate_keys_state=duplicate_keys)

    # Open and parse the JSON file, using the helper to track duplicates
    with open(json_file_path, 'r') as file:
        data = json.load(file, object_pairs_hook=hook)
        
    return _process_results(data, duplicate_keys)

def create_json_from_tsv(path_to_file, cell_type):

        # Validate evaluator input file exists
    if not os.path.exists(path_to_file):
        print(f"ERROR: Evaluator input file '{path_to_file}' not found.")
        raise FileNotFoundError(f"Evaluator input file not found: {path_to_file}")

    try:
        sequence_dataFrame = pd.read_csv(path_to_file, sep=' ', header=None)
        sequence_dataFrame.columns = ['variant', 'reference', 'alternate']

        variant_counts = sequence_dataFrame['variant'].value_counts()
        duplicates_with_counts = variant_counts[variant_counts > 1]

        print("Duplicated variants with counts:")
        print(duplicates_with_counts)   
    
        sequence_dataFrame_no_duplicates = sequence_dataFrame.drop_duplicates(subset='variant')
        print(sequence_dataFrame_no_duplicates.shape)
        # Create a dictionary, mapping each sequence name to its corresponding sequence
        variant_dict = {}

        for _, row in sequence_dataFrame_no_duplicates.iterrows():
            variant = row['variant']
            variant_dict[f"{variant}_reference"] = row['reference']
            variant_dict[f"{variant}_alternate"] = row['alternate']


        # Define the prediction tasks as a separate variable
        prediction_tasks = [
            {
                "name": "engritz_" + cell_type,
                "type": "expression",
                "cell_type": cell_type,
                "scale": "log",
                "species": "homo_sapiens"
            }
        ]
        
        # Build the JSON evaluator object
        evaluator_dict = {
            "request": "predict",
            "readout": "point",
            "prediction_tasks": prediction_tasks,
            "sequences": variant_dict
        }
        
        # Convert the dictionary to a JSON string with indentation for readability
        json_string = json.dumps(evaluator_dict)
        #Check for duplicate keys in the generated JSON string.
        # Use the helper function that accepts a JSON string.
        jsonResult_dict = check_duplicates_from_string(json_string)
        return jsonResult_dict
    except (json.JSONDecodeError, 
        DuplicateKeysError) as e:
        # Raise a general ValueError that the main script's handler
        # will catch and report cleanly
        raise ValueError(f"Input data is invalid.\nDetails: {e}") from e