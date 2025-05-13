import math
from utils import *
import time
import pandas as pd
from tqdm import trange
import os
from datetime import datetime
import csv
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-config', type=str, help="Path to the parameter_cfg.config", required=True)
parser.add_argument('-o', type=str, help="Path to the output file", required=True)
parser.add_argument('-n', type=int, default=3000, help="Sliding window size")
parser.add_argument('-k', type=float, default=0.3, help="The reference value for the allowed degree of deviation.")

args = parser.parse_args()


def calcu_win_depth(data, pgg, window_size):
    win_depth = []
    for i in range(len(data) - window_size + 1):
        x_bar = np.mean(data[i:i + window_size])
        p_bar = np.mean(pgg[i:i + window_size])
        s_bar = np.log(x_bar / p_bar)
        win_depth.append(s_bar)
    return win_depth


def calcu_ct(data, K):
    ct_up = np.zeros(len(data))
    ct_down = np.zeros(len(data))
    for i in range(len(data)):
        if i == 0:
            ct_up[i] = 0
            ct_down[i] = 0
        else:
            ct_up[i] = np.max([0, data[i] - K + ct_up[i - 1]])
            ct_down[i] = np.min([0, data[i] + K + ct_down[i - 1]])
    return ct_up, ct_down


def find_continuous_up_segments(ct, H_pos, min_length):
    up_segments = []
    current_segment = []

    for i, num in enumerate(ct):
        if num >= H_pos:
            current_segment.append((i, num))
        else:
            if len(current_segment) >= min_length:
                up_segments.append(current_segment)
            current_segment = []

    if len(current_segment) >= min_length:
        up_segments.append(current_segment)

    return up_segments


def find_continuous_down_segments(ct, H_neg, min_length):
    down_segments = []
    current_segment = []

    for i, num in enumerate(ct):
        if num <= H_neg:
            current_segment.append((i, num))
        else:
            if len(current_segment) >= min_length:
                down_segments.append(current_segment)
            current_segment = []

    if len(current_segment) >= min_length:
        down_segments.append(current_segment)

    return down_segments


def find_cand_dup_regs(result):
    cand_cnv_regs = []
    if result:
        for res in result:
            start = res[0][0]
            end = max(res, key=lambda item: item[1])[0]
            if end - start >= 1000:
                cand_cnv_regs.append([start, end])
        return cand_cnv_regs
    else:
        return []


def find_cand_del_regs(result):
    cand_cnv_regs = []
    if result:
        for res in result:
            start = res[0][0]
            end = min(res, key=lambda item: item[1])[0]
            if end - start >= 1000:
                cand_cnv_regs.append([start, end])
        return cand_cnv_regs
    else:
        return []


def calcu_logr(cand_cnv_regs, sample_depth, pgg_depth, chr_name, sample):
    det_results = []
    if cand_cnv_regs:
        for cand_cnv_reg in cand_cnv_regs:
            pgg_interval_depth = np.mean(pgg_depth[chr_name][cand_cnv_reg[0]: cand_cnv_reg[1]])
            sample_interval_depth = np.mean(sample_depth[chr_name][cand_cnv_reg[0]: cand_cnv_reg[1]])
            if sample_interval_depth != 0:
                logR = round(math.log(sample_interval_depth / pgg_interval_depth, 2), 5)
                if logR >= 0.3 or logR <= -0.3:
                    if logR > 0.3:
                        cnv_type = 'dup'
                    else:
                        cnv_type = 'del'
                    det_results.append([sample, chr_name, cand_cnv_reg[0], cand_cnv_reg[1], logR, cnv_type])
            else:
                logR = 'NA'
                cnv_type = 'del'
                det_results.append([sample, chr_name, cand_cnv_reg[0], cand_cnv_reg[1], logR, cnv_type])
    return det_results


def main():
    # .config file
    config_file_path = args.config
    paths = read_config(config_file_path)

    # The parameter of ZIP-Caller.
    # Please refer to the supplementary materials for setup details.
    slide_win = args.n
    K = args.k
    H_pos = np.log2(1.5) * slide_win
    H_neg = np.log2(0.5) * slide_win

    # Input: Sample files to be tested.
    test_file_list = paths['test_file_list']
    test_df = pd.read_csv(test_file_list, index_col=0)

    # Chromosome bed file
    chr_len_path = paths['chr_len_path']
    chr_len_list = read_chr_len_file(chr_len_path)

    # The path of baseline file
    baseline_save_path = paths['baseline_save_path']
    baseline_data = {}
    for chr_name, chr_len in chr_len_list:
        file_path = f'{baseline_save_path}/baseline_file_{chr_name}.npz'
        loaded_data = load_npz_file(file_path)
        baseline_data[chr_name] = loaded_data[chr_name]

    # Log
    log_filename = "log/ZIP-Caller_log.txt"
    os.makedirs(os.path.dirname(log_filename), exist_ok=True)

    if not os.path.exists(log_filename):
        with open(log_filename, 'w') as file:
            current_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            file.write(f"time: {current_time}\n")

    current_datetime = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

    # Output
    output_path = args.o
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    output_file = f'{output_path}/zipcaller_res_{current_datetime}.cnv'
    json_dir = 'data/nor/'

    # Start ZIP-Caller
    with open(output_file, 'w', newline='') as file:
        csv_writer = csv.writer(file, delimiter='\t')
        csv_writer.writerow(
            ['SampleID', 'Chromosome', 'Start', 'End', 'LogR_Ratio', 'CNV_Type'])
        for i in trange(test_df.shape[0]):
            time.sleep(0.01)
            mapping = f"sample_{i}"
            file_name = test_df.loc[test_df['mapping'] == mapping]['file_name'].values[0]
            print(f'process {file_name} ..................')
            sample = os.path.basename(file_name)
            filename_part_ = os.path.basename(file_name).replace('.bam', '.json')
            json_file = os.path.join(json_dir, filename_part_)
            print(json_file)
            if os.path.isfile(json_file):
                try:
                    standardized_depth = load_from_json(json_file)
                    for chr_name, chr_len in chr_len_list:
                        s = standardized_depth[chr_name]
                        b = baseline_data[chr_name]

                        cusum_statistic = calcu_win_depth(s, b, slide_win)
                        ct_up, ct_down = calcu_ct(cusum_statistic, K)
                        result_up = find_continuous_up_segments(ct_up, H_pos=H_pos, min_length=10000)
                        result_down = find_continuous_down_segments(ct_down, H_neg=H_neg, min_length=10000)

                        cand_dup_reg = find_cand_dup_regs(result_up)
                        cand_down_reg = find_cand_del_regs(result_down)

                        if cand_dup_reg:
                            dup_results = calcu_logr(cand_dup_reg, standardized_depth, baseline_data, chr_name, sample)
                            for result in dup_results:
                                csv_writer.writerow(result)
                        if cand_down_reg:
                            del_results = calcu_logr(cand_down_reg, standardized_depth, baseline_data, chr_name, sample)
                            for result in del_results:
                                csv_writer.writerow(result)
                except Exception as e:
                    current_time = datetime.now().strftime("%Y%m%d_%H%M%S")
                    error_message = f"{current_time}:Error occurred while reading {sample}: {e}"
                    print(error_message)
                    # Write the error message to the log file
                    with open(log_filename, "a") as log_file:
                        log_file.write(error_message + "\n")
                    continue


if __name__ == '__main__':
    st = time.time()
    main()
    et = time.time()
    rt = et - st
    print(f"Finish! runtime: {rt}sec")
