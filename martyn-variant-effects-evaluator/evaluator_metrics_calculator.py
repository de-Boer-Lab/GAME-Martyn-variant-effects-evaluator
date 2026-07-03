'''Calculate and save the final evaluation metrics.'''

# NOTE: Every evaluator will do this slightly differently depending on how the data is presented

import os
import sys
import json
import pandas as pd
import numpy as np
import itertools
from datetime import datetime, timezone
from scipy.stats import pearsonr
from config import EVALUATOR_NAME, EVALUATOR_INPUT_PATH


def calculate_and_save_metrics(cell_type, saved_predictions_path_current, output_dir):
    
    try:
        # Correlation calculation
        # NOTE: Every evaluator will do this slightly differently depending on how the data is presented  
        if os.path.exists(saved_predictions_path_current):
            print("Starting Evaluation Calculation and Saving as CSV for: " + cell_type)
            print(f"Using predictions from: {saved_predictions_path_current}")
        
        correlation_summary_filepath = f"correlation_summary_{EVALUATOR_NAME}_{cell_type}.csv"

        # Now load predictions
        try:
            with open(saved_predictions_path_current, 'r') as f:
                predictions_file_content = json.load(f)
        except Exception as e:
            print(f"FATAL: Could not load predictions from {saved_predictions_path_current}. {e}", file=sys.stderr)
            return

        # ADDITION: Construct file name after receiving predictor_name
        predictor_name_received = predictions_file_content.get("predictor_name", None)
        scale_prediction_actual = predictions_file_content['prediction_tasks'][0].get("scale_prediction_actual", None)

        predictor_name = predictor_name_received.replace(" ", "_").replace("/", "_")

        # Get UTC timestamp for predictor_nam
        timestamp = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S.%f")
        # Compute the full RETURN_FILE_PATH using the provided output directory

        print(f"Will save predictions to: {correlation_summary_filepath}")
        
        #Read in measured values
        measured_values_path = os.path.join(EVALUATOR_INPUT_PATH, cell_type, f"all_{cell_type}.tsv")

        print(f"Reading measured values from: {measured_values_path}")
        pearson_r = calculate_pearson_r(predictions_file_content, measured_values_path, scale_prediction_actual)

        description = f"Martyn Variant Effects ({cell_type})"
        #add code to create the output file

        prediction_task_data_onlyinfo = [{k: v for k, v in predictions_file_content["prediction_tasks"][0].items() if k != "predictions"}]

        evaluation_output = pd.DataFrame([{
            'Evaluator': EVALUATOR_NAME,
            'Description': description,
            'Predictor_name': predictor_name,
            'Time_stamp': timestamp,
            'Metric': 'pearson_r',
            'Value': str(pearson_r),
            'Prediction_task(s)_data': prediction_task_data_onlyinfo
        }])
        #print(evaluation_output)
        evaluation_output.to_csv(output_dir + '/' + correlation_summary_filepath , sep = "\t")
    
    except Exception as e:
        print(f"An unexpected error occurred during evaluation calculations: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()


def calculate_pearson_r(predictions_content, measured_values_path, scale_prediction_actual):
    """
    1- Calculates the log2FC between reference and alternate sequence predictions.
    2- Calculated the pearsonR between measured and predicted log2FC

    Args:
        predictions_json_path (str): Path to JSON file with predictions
    Returns:
        float: The Pearson correlation coefficient (r), or None if calculation isn't possible.
    """
   
    predictions_dict = predictions_content['prediction_tasks'][0]['predictions']
    if "error" in predictions_dict:
        print("No predictions were returned for this task -> Skipping evaluation calculation")
        return None

    # Create DataFrame from Predictions
    predictions_df = pd.DataFrame(list(predictions_dict.items()), columns=['id_column', 'Predicted_Value'])
    predictions_df['Predicted_Value'] = predictions_df['Predicted_Value'].apply(
        lambda x: x[0] if isinstance(x, list) and len(x) > 0 else x
        )
    print(predictions_df)
        
    #check here is there is NA is any of the prediction values
    na_rows = predictions_df[predictions_df['Predicted_Value'].isna()]
    if not na_rows.empty:
        print("NA values were found in the predictions, skipping evaluation")
        print(na_rows)
        return None
    
    # If model return predictions on linear scale log2fc = log2(alternate/reference)
    if scale_prediction_actual == 'linear':
        print("Warning: Predictor scale doesn't match the requested scale (log), be wary of the calculation")
        log2_fold_change_predicted = np.log2(predictions_df['Predicted_Value'].iloc[1::2].values/predictions_df['Predicted_Value'].iloc[0::2].values)

    # If model return predictions on log scale log2fc = log2(alternate) - log2(reference)
    if scale_prediction_actual == 'log':
        log2_fold_change_predicted =  predictions_df['Predicted_Value'].iloc[1::2].values - predictions_df['Predicted_Value'].iloc[0::2].values
    
    # Create a new DataFrame with results
    delta_df = pd.DataFrame({
        'ref_id': predictions_df['id_column'].iloc[0::2].values,
        'alt_id': predictions_df['id_column'].iloc[1::2].values,
        'log2_fold_change_predicted': log2_fold_change_predicted
    })

    delta_df['variant'] = delta_df['ref_id'].str.extract(r'^(.*)_(?:reference|alternate)$')
    print(delta_df)
    
    measured_values = pd.read_csv(measured_values_path, sep = '\t')
    print(measured_values)
 
    merged_data = measured_values.merge(delta_df, how='left', on='variant')
    print(merged_data.shape)
    print(merged_data[merged_data[['log2_fold_change', 'log2_fold_change_predicted']].isna().any(axis=1)])
    r, _ = pearsonr(merged_data['log2_fold_change'], merged_data['log2_fold_change_predicted'])
    print(f"Calculated Pearson r: {r}") 
    return r