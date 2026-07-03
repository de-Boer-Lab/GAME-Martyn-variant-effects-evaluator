'''Calculate and save the final evaluation metrics.'''

# NOTE: Every evaluator will do this slightly differently depending on how the data is presented

import os
import sys
import json
import numpy as np
import pandas as pd
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
        
        summary_filepath = os.path.join(output_dir, f"evaluation_summary_{EVALUATOR_NAME}.csv")

        # Now load predictions
        try:
            with open(saved_predictions_path_current, 'r') as f:
                predictions_file_content = json.load(f)
        except Exception as e:
            print(f"FATAL: Could not load predictions from {saved_predictions_path_current}. {e}", file=sys.stderr)
            return

        # CHANGE (v1): predictor_name fallback so a missing key doesn't AttributeError on .replace()
        predictor_name_received = predictions_file_content.get("predictor_name") or "UnknownPredictor"
        scale_prediction_actual = predictions_file_content['prediction_tasks'][0].get("scale_prediction_actual", None)

        predictor_name = predictor_name_received.replace(" ", "_").replace("/", "_")

        # Get UTC timestamp for predictor_nam
        timestamp = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S.%f")

        # print(f"Will save predictions to: {correlation_summary_filepath}")
        
        #Read in measured values
        measured_values_path = os.path.join(EVALUATOR_INPUT_PATH, cell_type, f"all_{cell_type}.tsv")
        print(f"Reading measured values from: {measured_values_path}")
        pearson_r = calculate_pearson_r(predictions_file_content, measured_values_path, scale_prediction_actual)

        description = f"Martyn Variant Effects ({cell_type})"
        #add code to create the output file

        prediction_task_data_onlyinfo = [{k: v for k, v in predictions_file_content["prediction_tasks"][0].items() if k != "predictions"}]
        
        # None (disqualified) -> "NaN" string; 0.0 (zero variance) stays 0.0.
        value_str = "NaN" if pearson_r is None else str(pearson_r)

        evaluation_output = pd.DataFrame([{
            'evaluator_name': EVALUATOR_NAME,
            'description': description,
            'predictor_name': predictor_name,
            'time_stamp': timestamp,
            'metric': 'pearson_r',
            'value': value_str,
            'prediction_task(s)_data': prediction_task_data_onlyinfo
        }])
        
        # Append to the single summary file (header only when new), no index.
        file_exists = os.path.isfile(summary_filepath)
        evaluation_output.to_csv(summary_filepath, mode='a', sep='\t',
                                 header=(not file_exists), index=False)
        if file_exists:
            print(f"Appended metrics to {summary_filepath}")
        else:
            print(f"Created new metrics file {summary_filepath}")
    
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
        measured_values_path (str): Path to TSV file with measured values
        scale_prediction_actual (str): The scale of the predictions, either 'log' or 'linear'
        
    Returns:
        float: The Pearson correlation coefficient (r), or None if calculation isn't possible.
        0.0  : If predicted or measured values have zero variance ("ran but useless").
        None : If disqualified (no/NA predictions, wrong scale, nothing to correlate).
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
    
    predictions_df.columns = ['VariantID_h19', 'Predicted_Value']
    print(predictions_df)
    #check here is there is NA is any of the prediction values
    na_rows = predictions_df[predictions_df['Predicted_Value'].isna()]
    if not na_rows.empty:
        print("NA values were found in the predictions, skipping evaluation")
        print(na_rows)
        return None

    measured_values = pd.read_csv(measured_values_path, sep = '\t')
    print(measured_values)
    
    merged_data = measured_values.merge(predictions_df, how='inner', on='VariantID_h19')
    print(merged_data.shape)

    # Get PPIF value
    ppif_pred = predictions_df.loc[predictions_df['VariantID_h19'] == "PPIF", 'Predicted_Value']
    ppif_pred = ppif_pred.iloc[0] if not ppif_pred.empty else None

    # Get IL2RA value
    il2ra_pred = predictions_df.loc[predictions_df['VariantID_h19'] == "IL2RA", 'Predicted_Value']
    il2ra_pred = il2ra_pred.iloc[0] if not il2ra_pred.empty else None

    print("PPIF Reference Prediction:", ppif_pred)
    print("IL2RA Reference Prediction:", il2ra_pred)

    # Compute predicted log2 fold change per variant, relative to its gene's reference prediction
    # log case:    predictions are already log2-scaled -> log2FC = pred_alt - pred_ref
    # linear case: predictions are linear              -> log2FC = log2(pred_alt / pred_ref)
    # Both recover the same quantity (and the same Pearson r, which is invariant to log base),
    # so a predictor that returns linear is still evaluated correctly rather than rejected
    if scale_prediction_actual == 'log':
        for idx, row in merged_data.iterrows():
            gene = row['gene_symbol']
            if gene == "PPIF" and ppif_pred is not None:
                merged_data.loc[idx, "Predicted_Value_log2FC"] = row["Predicted_Value"] - ppif_pred
            if gene == "IL2RA" and il2ra_pred is not None:
                merged_data.loc[idx, "Predicted_Value_log2FC"] = row["Predicted_Value"] - il2ra_pred

    elif scale_prediction_actual == 'linear':
        # CHANGE (v1): was `row["Predicted_Value"] / 1` -- a no-op placeholder that left the value
        # on a linear scale and silently mixed scales into the correlation. log2(alt/ref) is the
        # correct conversion (matches the original code comment). Guarded for positivity: log2 of
        # a non-positive ratio is undefined, so those rows are left NaN and dropped below
        print("Warning: Predictor returned linear (requested log). Converting via log2(alt/ref).")
        for idx, row in merged_data.iterrows():
            gene = row['gene_symbol']
            alt_pred = row["Predicted_Value"]
            if gene == "PPIF" and ppif_pred is not None and ppif_pred > 0 and alt_pred > 0:
                merged_data.loc[idx, "Predicted_Value_log2FC"] = np.log2(alt_pred / ppif_pred)
            if gene == "IL2RA" and il2ra_pred is not None and il2ra_pred > 0 and alt_pred > 0:
                merged_data.loc[idx, "Predicted_Value_log2FC"] = np.log2(alt_pred / il2ra_pred)


    print(merged_data)
    
    if "Predicted_Value_log2FC" not in merged_data.columns:
        print("No log2FC values were computed (no PPIF/IL2RA variants matched). Skipping.")
        return None
    
    merged_data['log2_fold_change'] = pd.to_numeric(merged_data['log2_fold_change'], errors='coerce')
    merged_data['Predicted_Value_log2FC'] = pd.to_numeric(merged_data['Predicted_Value_log2FC'], errors='coerce')
    before = len(merged_data)
    merged_data = merged_data.dropna(subset=['log2_fold_change', 'Predicted_Value_log2FC'])
    dropped = before - len(merged_data)
    if dropped:
        print(f"Dropped {dropped} row(s) with non-numeric/missing values before correlation.")

    if len(merged_data) < 2:
        print("Fewer than 2 paired values after merge -> cannot compute correlation. Skipping.")
        return None

    # CHANGE (v1): zero-variance handling (pearsonr is undefined for a constant input).
    if merged_data['log2_fold_change'].std() == 0 or merged_data['Predicted_Value_log2FC'].std() == 0:
        print("Zero variance in measured or predicted log2FC -> assigning pearson_r = 0.0")
        return 0.0

    try:
        r, _ = pearsonr(merged_data['log2_fold_change'], merged_data['Predicted_Value_log2FC'])
        pearson_r = 0.0 if np.isnan(r) else float(r)   # CHANGE (v1): NaN result -> 0.0
        print(f"Calculated Pearson r: {pearson_r}")
        return pearson_r
    except Exception as e:
        print(f"An error occurred during the correlation calculation: {e}", file=sys.stderr)
        return None

