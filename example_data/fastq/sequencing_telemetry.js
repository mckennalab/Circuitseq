[
    {
        "aggregation": "segment", 
        "analysis_id": "e2692740-f853-4dfa-99c2-3309f56cdfdb", 
        "basecall_1d": {
            "exit_status_dist": {
                "pass": 12
            }, 
            "qscore_dist_temp": [
                {
                    "count": 1, 
                    "mean_qscore": 10.0
                }, 
                {
                    "count": 1, 
                    "mean_qscore": 11.0
                }, 
                {
                    "count": 4, 
                    "mean_qscore": 12.5
                }, 
                {
                    "count": 2, 
                    "mean_qscore": 13.5
                }, 
                {
                    "count": 2, 
                    "mean_qscore": 14.0
                }, 
                {
                    "count": 2, 
                    "mean_qscore": 14.5
                }
            ], 
            "qscore_sum_temp": {
                "count": 12, 
                "mean": 13.230557441711426, 
                "sum": 158.76669311523438
            }, 
            "read_len_events_sum_temp": 282516, 
            "seq_len_bases_dist_temp": [
                {
                    "count": 12, 
                    "length": 0.0
                }
            ], 
            "seq_len_bases_sum_temp": 12, 
            "seq_len_events_dist_temp": [
                {
                    "count": 1, 
                    "length": 16000.0
                }, 
                {
                    "count": 1, 
                    "length": 18000.0
                }, 
                {
                    "count": 1, 
                    "length": 19000.0
                }, 
                {
                    "count": 1, 
                    "length": 21000.0
                }, 
                {
                    "count": 3, 
                    "length": 23000.0
                }, 
                {
                    "count": 3, 
                    "length": 25000.0
                }, 
                {
                    "count": 1, 
                    "length": 27000.0
                }, 
                {
                    "count": 1, 
                    "length": 31000.0
                }
            ], 
            "speed_bases_per_second_dist_temp": [
                {
                    "count": 12, 
                    "speed": 1.0
                }
            ], 
            "strand_median_pa": {
                "count": 12, 
                "mean": 77.33541107177734, 
                "sum": 928.02490234375
            }, 
            "strand_sd_pa": {
                "count": 12, 
                "mean": 9.652189254760742, 
                "sum": 115.8262710571289
            }
        }, 
        "channel_count": 12, 
        "context_tags": {
            "barcoding_enabled": "0", 
            "experiment_duration_set": "1680", 
            "experiment_type": "genomic_dna", 
            "local_basecalling": "0", 
            "package": "bream4", 
            "package_version": "6.1.10", 
            "sample_frequency": "4000", 
            "sequencing_kit": "sqk-lsk110"
        }, 
        "latest_run_time": 3471.63720703125, 
        "levels_sums": {
            "count": 12, 
            "mean": 204.1630859375, 
            "open_pore_level_sum": 2449.95703125
        }, 
        "opts": {
            "adapter_pt_range_scale": "5.200000", 
            "additional_lamp_context_bases": "2", 
            "align_ref": "", 
            "align_type": "auto", 
            "allow_inferior_barcodes": "0", 
            "as_cpu_threads_per_scaler": "2", 
            "as_gpu_runners_per_device": "2", 
            "as_model_file": "", 
            "as_num_scalers": "4", 
            "as_reads_per_runner": "32", 
            "bam_methylation_threshold": "5.000000", 
            "bam_out": "0", 
            "barcode_kits": "", 
            "barcode_nested_output_folder": "0", 
            "beam_cut": "100.000000", 
            "beam_width": "32", 
            "bed_file": "", 
            "builtin_scripts": "1", 
            "calib_detect": "0", 
            "calib_max_sequence_length": "3800", 
            "calib_min_coverage": "0.600000", 
            "calib_min_sequence_length": "3000", 
            "calib_reference": "lambda_3.6kb.fasta", 
            "chunk_size": "2000", 
            "chunks_per_caller": "10000", 
            "chunks_per_runner": "208", 
            "client_id": "-1", 
            "compress_fastq": "1", 
            "cpu_threads_per_caller": "4", 
            "detect_adapter": "0", 
            "detect_barcodes": "0", 
            "detect_mid_strand_adapter": "0", 
            "detect_mid_strand_barcodes": "0", 
            "detect_primer": "0", 
            "device": "cuda:0", 
            "disable_pings": "0", 
            "disable_qscore_filtering": "0", 
            "dmean_threshold": "1.000000", 
            "dmean_win_size": "2", 
            "do_read_splitting": "0", 
            "duplex_window_size_max": "1000", 
            "duplex_window_size_min": "200", 
            "end_gap1": "40", 
            "end_gap2": "40", 
            "extend_gap1": "40", 
            "extend_gap2": "160", 
            "fast5_out": "0", 
            "flowcell": "", 
            "front_window_size": "150", 
            "gpu_runners_per_device": "12", 
            "high_priority_threshold": "10", 
            "index": "0", 
            "input_file_list": "", 
            "int8_mode": "0", 
            "jump_threshold": "1.000000", 
            "kernel_path": "", 
            "kit": "", 
            "lamp_kit": "", 
            "log_speed_frequency": "0", 
            "max_queued_reads": "2000", 
            "max_read_split_depth": "2", 
            "max_search_len": "1000", 
            "medium_priority_threshold": "4", 
            "min_length_lamp_context": "30", 
            "min_length_lamp_target": "70", 
            "min_qscore": "10.000000", 
            "min_score_adapter": "60.000000", 
            "min_score_adapter_mid": "50.000000", 
            "min_score_barcode_front": "60.000000", 
            "min_score_barcode_mask": "40.000000", 
            "min_score_barcode_mid": "50.000000", 
            "min_score_barcode_rear": "60.000000", 
            "min_score_lamp": "80.000000", 
            "min_score_lamp_mask": "50.000000", 
            "min_score_lamp_target": "50.000000", 
            "min_score_primer": "60.000000", 
            "min_score_read_splitting": "70.000000", 
            "model_file": "template_r9.4.1_450bps_sup.jsn", 
            "moves_out": "0", 
            "nested_output_folder": "0", 
            "noisiest_section_scaling_max_size": "8000", 
            "num_alignment_threads": "4", 
            "num_barcode_threads": "4", 
            "num_barcoding_buffers": "24", 
            "num_base_mod_threads": "2", 
            "num_callers": "4", 
            "num_extra_bases_trim": "0", 
            "num_mid_barcoding_buffers": "96", 
            "num_read_splitting_buffers": "16", 
            "num_read_splitting_threads": "4", 
            "num_reads_per_barcoding_buffer": "4", 
            "open_gap1": "40", 
            "open_gap2": "160", 
            "overlap": "100", 
            "override_scaling": "0", 
            "ping_segment_duration": "60", 
            "ping_url": "https://ping.oxfordnanoportal.com/basecall", 
            "post_out": "0", 
            "print_workflows": "0", 
            "progress_stats_frequency": "-1.000000", 
            "pt_median_offset": "2.500000", 
            "pt_minimum_read_start_index": "30", 
            "pt_required_adapter_drop": "30.000000", 
            "pt_scaling": "0", 
            "qscore_offset": "0.349800", 
            "qscore_scale": "0.972200", 
            "quiet": "0", 
            "read_batch_size": "4000", 
            "read_id_list": "", 
            "read_splitting_arrangement_files": "", 
            "read_splitting_score_matrix_filename": "", 
            "rear_window_size": "150", 
            "records_per_fastq": "0", 
            "recursive": "0", 
            "remora_models": "", 
            "require_barcodes_both_ends": "0", 
            "resume": "0", 
            "reverse_sequence": "0", 
            "sample_sheet": "", 
            "scaling_mad": "1.000000", 
            "scaling_med": "0.000000", 
            "score_matrix_filename": "5x5_mismatch_matrix.txt", 
            "start_gap1": "40", 
            "start_gap2": "40", 
            "stay_penalty": "1.000000", 
            "temp_bias": "1.000000", 
            "temp_weight": "1.000000", 
            "trace_categories_logs": "", 
            "trace_domains_config": "", 
            "trim_adapters": "0", 
            "trim_barcodes": "0", 
            "trim_min_events": "3", 
            "trim_primers": "0", 
            "trim_strategy": "dna", 
            "trim_threshold": "2.500000", 
            "u_substitution": "0", 
            "verbose_logs": "0"
        }, 
        "read_count": 12, 
        "reads_per_channel_dist": [
            {
                "channel": 9, 
                "count": 1
            }, 
            {
                "channel": 16, 
                "count": 1
            }, 
            {
                "channel": 27, 
                "count": 1
            }, 
            {
                "channel": 28, 
                "count": 1
            }, 
            {
                "channel": 31, 
                "count": 1
            }, 
            {
                "channel": 49, 
                "count": 1
            }, 
            {
                "channel": 54, 
                "count": 1
            }, 
            {
                "channel": 67, 
                "count": 1
            }, 
            {
                "channel": 71, 
                "count": 1
            }, 
            {
                "channel": 99, 
                "count": 1
            }, 
            {
                "channel": 100, 
                "count": 1
            }, 
            {
                "channel": 115, 
                "count": 1
            }
        ], 
        "run_id": "6a94eb36563a9975b38e81c41534f6590c24867c", 
        "segment_duration": 60, 
        "segment_number": 1, 
        "segment_type": "guppy-acquisition", 
        "software": {
            "analysis": "1d_basecalling", 
            "name": "guppy-basecalling", 
            "version": "6.1.2+e0556ff"
        }, 
        "tracking_id": {
            "asic_id": "3506339664", 
            "asic_id_eeprom": "5924305", 
            "asic_temp": "35.086052", 
            "asic_version": "IA02D", 
            "auto_update": "0", 
            "auto_update_source": "https://mirror.oxfordnanoportal.com/software/MinKNOW/", 
            "bream_is_standard": "0", 
            "configuration_version": "4.2.11", 
            "device_id": "MN34875", 
            "device_type": "minion", 
            "distribution_status": "stable", 
            "distribution_version": "21.02.1", 
            "exp_script_name": "sequencing_MIN106_DNA:FLO-FLG001:SQK-LSK110", 
            "exp_script_purpose": "sequencing_run", 
            "exp_start_time": "2021-08-03T19:16:19Z", 
            "flongle_adapter_id": "FA-00993", 
            "flow_cell_id": "AHC239", 
            "flow_cell_product_code": "FLO-FLG001", 
            "guppy_version": "4.3.4+ecb2805", 
            "heatsink_temp": "35.003906", 
            "hostname": "mckennalab", 
            "installation_type": "nc", 
            "local_firmware_file": "1", 
            "msg_id": "743d916d-1b36-46b4-94dc-612918f0ab87", 
            "operating_system": "ubuntu 20.04", 
            "protocol_group_id": "FEX_120", 
            "protocol_run_id": "bdcd3859-b6c1-400c-b4e6-2bc287733f75", 
            "protocol_start_time": "", 
            "protocols_version": "6.1.10", 
            "run_id": "6a94eb36563a9975b38e81c41534f6590c24867c", 
            "sample_id": "96_v2", 
            "time_stamp": "2022-06-11T02:04:21Z", 
            "usb_config": "MinION_fx3_1.1.1_ONT#MinION_fpga_1.1.0#bulk#Auto", 
            "version": "4.2.5"
        }
    }, 
    {
        "aggregation": "segment", 
        "analysis_id": "e2692740-f853-4dfa-99c2-3309f56cdfdb", 
        "basecall_1d": {
            "exit_status_dist": {
                "pass": 8
            }, 
            "qscore_dist_temp": [
                {
                    "count": 1, 
                    "mean_qscore": 10.5
                }, 
                {
                    "count": 1, 
                    "mean_qscore": 11.0
                }, 
                {
                    "count": 1, 
                    "mean_qscore": 12.0
                }, 
                {
                    "count": 2, 
                    "mean_qscore": 12.5
                }, 
                {
                    "count": 1, 
                    "mean_qscore": 13.0
                }, 
                {
                    "count": 1, 
                    "mean_qscore": 13.5
                }, 
                {
                    "count": 1, 
                    "mean_qscore": 15.0
                }
            ], 
            "qscore_sum_temp": {
                "count": 8, 
                "mean": 12.767404556274414, 
                "sum": 102.13923645019531
            }, 
            "read_len_events_sum_temp": 212069, 
            "seq_len_bases_dist_temp": [
                {
                    "count": 8, 
                    "length": 0.0
                }
            ], 
            "seq_len_bases_sum_temp": 8, 
            "seq_len_events_dist_temp": [
                {
                    "count": 1, 
                    "length": 20000.0
                }, 
                {
                    "count": 1, 
                    "length": 23000.0
                }, 
                {
                    "count": 1, 
                    "length": 24000.0
                }, 
                {
                    "count": 1, 
                    "length": 25000.0
                }, 
                {
                    "count": 2, 
                    "length": 28000.0
                }, 
                {
                    "count": 1, 
                    "length": 30000.0
                }, 
                {
                    "count": 1, 
                    "length": 31000.0
                }
            ], 
            "speed_bases_per_second_dist_temp": [
                {
                    "count": 8, 
                    "speed": 1.0
                }
            ], 
            "strand_median_pa": {
                "count": 8, 
                "mean": 73.25349426269531, 
                "sum": 586.0279541015625
            }, 
            "strand_sd_pa": {
                "count": 8, 
                "mean": 9.305890083312988, 
                "sum": 74.4471206665039
            }
        }, 
        "channel_count": 8, 
        "context_tags": {
            "barcoding_enabled": "0", 
            "experiment_duration_set": "1680", 
            "experiment_type": "genomic_dna", 
            "local_basecalling": "0", 
            "package": "bream4", 
            "package_version": "6.1.10", 
            "sample_frequency": "4000", 
            "sequencing_kit": "sqk-lsk110"
        }, 
        "latest_run_time": 6941.009765625, 
        "levels_sums": {
            "count": 8, 
            "mean": 198.3103485107422, 
            "open_pore_level_sum": 1586.4827880859375
        }, 
        "opts": {
            "adapter_pt_range_scale": "5.200000", 
            "additional_lamp_context_bases": "2", 
            "align_ref": "", 
            "align_type": "auto", 
            "allow_inferior_barcodes": "0", 
            "as_cpu_threads_per_scaler": "2", 
            "as_gpu_runners_per_device": "2", 
            "as_model_file": "", 
            "as_num_scalers": "4", 
            "as_reads_per_runner": "32", 
            "bam_methylation_threshold": "5.000000", 
            "bam_out": "0", 
            "barcode_kits": "", 
            "barcode_nested_output_folder": "0", 
            "beam_cut": "100.000000", 
            "beam_width": "32", 
            "bed_file": "", 
            "builtin_scripts": "1", 
            "calib_detect": "0", 
            "calib_max_sequence_length": "3800", 
            "calib_min_coverage": "0.600000", 
            "calib_min_sequence_length": "3000", 
            "calib_reference": "lambda_3.6kb.fasta", 
            "chunk_size": "2000", 
            "chunks_per_caller": "10000", 
            "chunks_per_runner": "208", 
            "client_id": "-1", 
            "compress_fastq": "1", 
            "cpu_threads_per_caller": "4", 
            "detect_adapter": "0", 
            "detect_barcodes": "0", 
            "detect_mid_strand_adapter": "0", 
            "detect_mid_strand_barcodes": "0", 
            "detect_primer": "0", 
            "device": "cuda:0", 
            "disable_pings": "0", 
            "disable_qscore_filtering": "0", 
            "dmean_threshold": "1.000000", 
            "dmean_win_size": "2", 
            "do_read_splitting": "0", 
            "duplex_window_size_max": "1000", 
            "duplex_window_size_min": "200", 
            "end_gap1": "40", 
            "end_gap2": "40", 
            "extend_gap1": "40", 
            "extend_gap2": "160", 
            "fast5_out": "0", 
            "flowcell": "", 
            "front_window_size": "150", 
            "gpu_runners_per_device": "12", 
            "high_priority_threshold": "10", 
            "index": "0", 
            "input_file_list": "", 
            "int8_mode": "0", 
            "jump_threshold": "1.000000", 
            "kernel_path": "", 
            "kit": "", 
            "lamp_kit": "", 
            "log_speed_frequency": "0", 
            "max_queued_reads": "2000", 
            "max_read_split_depth": "2", 
            "max_search_len": "1000", 
            "medium_priority_threshold": "4", 
            "min_length_lamp_context": "30", 
            "min_length_lamp_target": "70", 
            "min_qscore": "10.000000", 
            "min_score_adapter": "60.000000", 
            "min_score_adapter_mid": "50.000000", 
            "min_score_barcode_front": "60.000000", 
            "min_score_barcode_mask": "40.000000", 
            "min_score_barcode_mid": "50.000000", 
            "min_score_barcode_rear": "60.000000", 
            "min_score_lamp": "80.000000", 
            "min_score_lamp_mask": "50.000000", 
            "min_score_lamp_target": "50.000000", 
            "min_score_primer": "60.000000", 
            "min_score_read_splitting": "70.000000", 
            "model_file": "template_r9.4.1_450bps_sup.jsn", 
            "moves_out": "0", 
            "nested_output_folder": "0", 
            "noisiest_section_scaling_max_size": "8000", 
            "num_alignment_threads": "4", 
            "num_barcode_threads": "4", 
            "num_barcoding_buffers": "24", 
            "num_base_mod_threads": "2", 
            "num_callers": "4", 
            "num_extra_bases_trim": "0", 
            "num_mid_barcoding_buffers": "96", 
            "num_read_splitting_buffers": "16", 
            "num_read_splitting_threads": "4", 
            "num_reads_per_barcoding_buffer": "4", 
            "open_gap1": "40", 
            "open_gap2": "160", 
            "overlap": "100", 
            "override_scaling": "0", 
            "ping_segment_duration": "60", 
            "ping_url": "https://ping.oxfordnanoportal.com/basecall", 
            "post_out": "0", 
            "print_workflows": "0", 
            "progress_stats_frequency": "-1.000000", 
            "pt_median_offset": "2.500000", 
            "pt_minimum_read_start_index": "30", 
            "pt_required_adapter_drop": "30.000000", 
            "pt_scaling": "0", 
            "qscore_offset": "0.349800", 
            "qscore_scale": "0.972200", 
            "quiet": "0", 
            "read_batch_size": "4000", 
            "read_id_list": "", 
            "read_splitting_arrangement_files": "", 
            "read_splitting_score_matrix_filename": "", 
            "rear_window_size": "150", 
            "records_per_fastq": "0", 
            "recursive": "0", 
            "remora_models": "", 
            "require_barcodes_both_ends": "0", 
            "resume": "0", 
            "reverse_sequence": "0", 
            "sample_sheet": "", 
            "scaling_mad": "1.000000", 
            "scaling_med": "0.000000", 
            "score_matrix_filename": "5x5_mismatch_matrix.txt", 
            "start_gap1": "40", 
            "start_gap2": "40", 
            "stay_penalty": "1.000000", 
            "temp_bias": "1.000000", 
            "temp_weight": "1.000000", 
            "trace_categories_logs": "", 
            "trace_domains_config": "", 
            "trim_adapters": "0", 
            "trim_barcodes": "0", 
            "trim_min_events": "3", 
            "trim_primers": "0", 
            "trim_strategy": "dna", 
            "trim_threshold": "2.500000", 
            "u_substitution": "0", 
            "verbose_logs": "0"
        }, 
        "read_count": 8, 
        "reads_per_channel_dist": [
            {
                "channel": 9, 
                "count": 1
            }, 
            {
                "channel": 18, 
                "count": 1
            }, 
            {
                "channel": 26, 
                "count": 1
            }, 
            {
                "channel": 31, 
                "count": 1
            }, 
            {
                "channel": 32, 
                "count": 1
            }, 
            {
                "channel": 34, 
                "count": 1
            }, 
            {
                "channel": 82, 
                "count": 1
            }, 
            {
                "channel": 89, 
                "count": 1
            }
        ], 
        "run_id": "6a94eb36563a9975b38e81c41534f6590c24867c", 
        "segment_duration": 60, 
        "segment_number": 2, 
        "segment_type": "guppy-acquisition", 
        "software": {
            "analysis": "1d_basecalling", 
            "name": "guppy-basecalling", 
            "version": "6.1.2+e0556ff"
        }, 
        "tracking_id": {
            "asic_id": "3506339664", 
            "asic_id_eeprom": "5924305", 
            "asic_temp": "35.086052", 
            "asic_version": "IA02D", 
            "auto_update": "0", 
            "auto_update_source": "https://mirror.oxfordnanoportal.com/software/MinKNOW/", 
            "bream_is_standard": "0", 
            "configuration_version": "4.2.11", 
            "device_id": "MN34875", 
            "device_type": "minion", 
            "distribution_status": "stable", 
            "distribution_version": "21.02.1", 
            "exp_script_name": "sequencing_MIN106_DNA:FLO-FLG001:SQK-LSK110", 
            "exp_script_purpose": "sequencing_run", 
            "exp_start_time": "2021-08-03T19:16:19Z", 
            "flongle_adapter_id": "FA-00993", 
            "flow_cell_id": "AHC239", 
            "flow_cell_product_code": "FLO-FLG001", 
            "guppy_version": "4.3.4+ecb2805", 
            "heatsink_temp": "35.003906", 
            "hostname": "mckennalab", 
            "installation_type": "nc", 
            "local_firmware_file": "1", 
            "msg_id": "482cee99-28d3-4465-809a-0996568097bf", 
            "operating_system": "ubuntu 20.04", 
            "protocol_group_id": "FEX_120", 
            "protocol_run_id": "bdcd3859-b6c1-400c-b4e6-2bc287733f75", 
            "protocol_start_time": "", 
            "protocols_version": "6.1.10", 
            "run_id": "6a94eb36563a9975b38e81c41534f6590c24867c", 
            "sample_id": "96_v2", 
            "time_stamp": "2022-06-11T02:04:21Z", 
            "usb_config": "MinION_fx3_1.1.1_ONT#MinION_fpga_1.1.0#bulk#Auto", 
            "version": "4.2.5"
        }
    }, 
    {
        "aggregation": "segment", 
        "analysis_id": "e2692740-f853-4dfa-99c2-3309f56cdfdb", 
        "basecall_1d": {
            "exit_status_dist": {
                "pass": 7
            }, 
            "qscore_dist_temp": [
                {
                    "count": 1, 
                    "mean_qscore": 10.5
                }, 
                {
                    "count": 1, 
                    "mean_qscore": 11.0
                }, 
                {
                    "count": 1, 
                    "mean_qscore": 11.5
                }, 
                {
                    "count": 2, 
                    "mean_qscore": 12.5
                }, 
                {
                    "count": 1, 
                    "mean_qscore": 14.0
                }, 
                {
                    "count": 1, 
                    "mean_qscore": 14.5
                }
            ], 
            "qscore_sum_temp": {
                "count": 7, 
                "mean": 12.586404800415039, 
                "sum": 88.1048355102539
            }, 
            "read_len_events_sum_temp": 168095, 
            "seq_len_bases_dist_temp": [
                {
                    "count": 7, 
                    "length": 0.0
                }
            ], 
            "seq_len_bases_sum_temp": 7, 
            "seq_len_events_dist_temp": [
                {
                    "count": 1, 
                    "length": 16000.0
                }, 
                {
                    "count": 2, 
                    "length": 23000.0
                }, 
                {
                    "count": 1, 
                    "length": 24000.0
                }, 
                {
                    "count": 2, 
                    "length": 26000.0
                }, 
                {
                    "count": 1, 
                    "length": 27000.0
                }
            ], 
            "speed_bases_per_second_dist_temp": [
                {
                    "count": 7, 
                    "speed": 1.0
                }
            ], 
            "strand_median_pa": {
                "count": 7, 
                "mean": 70.45572662353516, 
                "sum": 493.1900634765625
            }, 
            "strand_sd_pa": {
                "count": 7, 
                "mean": 9.094321250915527, 
                "sum": 63.660247802734375
            }
        }, 
        "channel_count": 7, 
        "context_tags": {
            "barcoding_enabled": "0", 
            "experiment_duration_set": "1680", 
            "experiment_type": "genomic_dna", 
            "local_basecalling": "0", 
            "package": "bream4", 
            "package_version": "6.1.10", 
            "sample_frequency": "4000", 
            "sequencing_kit": "sqk-lsk110"
        }, 
        "latest_run_time": 9779.5234375, 
        "levels_sums": {
            "count": 7, 
            "mean": 191.9488067626953, 
            "open_pore_level_sum": 1343.6416015625
        }, 
        "opts": {
            "adapter_pt_range_scale": "5.200000", 
            "additional_lamp_context_bases": "2", 
            "align_ref": "", 
            "align_type": "auto", 
            "allow_inferior_barcodes": "0", 
            "as_cpu_threads_per_scaler": "2", 
            "as_gpu_runners_per_device": "2", 
            "as_model_file": "", 
            "as_num_scalers": "4", 
            "as_reads_per_runner": "32", 
            "bam_methylation_threshold": "5.000000", 
            "bam_out": "0", 
            "barcode_kits": "", 
            "barcode_nested_output_folder": "0", 
            "beam_cut": "100.000000", 
            "beam_width": "32", 
            "bed_file": "", 
            "builtin_scripts": "1", 
            "calib_detect": "0", 
            "calib_max_sequence_length": "3800", 
            "calib_min_coverage": "0.600000", 
            "calib_min_sequence_length": "3000", 
            "calib_reference": "lambda_3.6kb.fasta", 
            "chunk_size": "2000", 
            "chunks_per_caller": "10000", 
            "chunks_per_runner": "208", 
            "client_id": "-1", 
            "compress_fastq": "1", 
            "cpu_threads_per_caller": "4", 
            "detect_adapter": "0", 
            "detect_barcodes": "0", 
            "detect_mid_strand_adapter": "0", 
            "detect_mid_strand_barcodes": "0", 
            "detect_primer": "0", 
            "device": "cuda:0", 
            "disable_pings": "0", 
            "disable_qscore_filtering": "0", 
            "dmean_threshold": "1.000000", 
            "dmean_win_size": "2", 
            "do_read_splitting": "0", 
            "duplex_window_size_max": "1000", 
            "duplex_window_size_min": "200", 
            "end_gap1": "40", 
            "end_gap2": "40", 
            "extend_gap1": "40", 
            "extend_gap2": "160", 
            "fast5_out": "0", 
            "flowcell": "", 
            "front_window_size": "150", 
            "gpu_runners_per_device": "12", 
            "high_priority_threshold": "10", 
            "index": "0", 
            "input_file_list": "", 
            "int8_mode": "0", 
            "jump_threshold": "1.000000", 
            "kernel_path": "", 
            "kit": "", 
            "lamp_kit": "", 
            "log_speed_frequency": "0", 
            "max_queued_reads": "2000", 
            "max_read_split_depth": "2", 
            "max_search_len": "1000", 
            "medium_priority_threshold": "4", 
            "min_length_lamp_context": "30", 
            "min_length_lamp_target": "70", 
            "min_qscore": "10.000000", 
            "min_score_adapter": "60.000000", 
            "min_score_adapter_mid": "50.000000", 
            "min_score_barcode_front": "60.000000", 
            "min_score_barcode_mask": "40.000000", 
            "min_score_barcode_mid": "50.000000", 
            "min_score_barcode_rear": "60.000000", 
            "min_score_lamp": "80.000000", 
            "min_score_lamp_mask": "50.000000", 
            "min_score_lamp_target": "50.000000", 
            "min_score_primer": "60.000000", 
            "min_score_read_splitting": "70.000000", 
            "model_file": "template_r9.4.1_450bps_sup.jsn", 
            "moves_out": "0", 
            "nested_output_folder": "0", 
            "noisiest_section_scaling_max_size": "8000", 
            "num_alignment_threads": "4", 
            "num_barcode_threads": "4", 
            "num_barcoding_buffers": "24", 
            "num_base_mod_threads": "2", 
            "num_callers": "4", 
            "num_extra_bases_trim": "0", 
            "num_mid_barcoding_buffers": "96", 
            "num_read_splitting_buffers": "16", 
            "num_read_splitting_threads": "4", 
            "num_reads_per_barcoding_buffer": "4", 
            "open_gap1": "40", 
            "open_gap2": "160", 
            "overlap": "100", 
            "override_scaling": "0", 
            "ping_segment_duration": "60", 
            "ping_url": "https://ping.oxfordnanoportal.com/basecall", 
            "post_out": "0", 
            "print_workflows": "0", 
            "progress_stats_frequency": "-1.000000", 
            "pt_median_offset": "2.500000", 
            "pt_minimum_read_start_index": "30", 
            "pt_required_adapter_drop": "30.000000", 
            "pt_scaling": "0", 
            "qscore_offset": "0.349800", 
            "qscore_scale": "0.972200", 
            "quiet": "0", 
            "read_batch_size": "4000", 
            "read_id_list": "", 
            "read_splitting_arrangement_files": "", 
            "read_splitting_score_matrix_filename": "", 
            "rear_window_size": "150", 
            "records_per_fastq": "0", 
            "recursive": "0", 
            "remora_models": "", 
            "require_barcodes_both_ends": "0", 
            "resume": "0", 
            "reverse_sequence": "0", 
            "sample_sheet": "", 
            "scaling_mad": "1.000000", 
            "scaling_med": "0.000000", 
            "score_matrix_filename": "5x5_mismatch_matrix.txt", 
            "start_gap1": "40", 
            "start_gap2": "40", 
            "stay_penalty": "1.000000", 
            "temp_bias": "1.000000", 
            "temp_weight": "1.000000", 
            "trace_categories_logs": "", 
            "trace_domains_config": "", 
            "trim_adapters": "0", 
            "trim_barcodes": "0", 
            "trim_min_events": "3", 
            "trim_primers": "0", 
            "trim_strategy": "dna", 
            "trim_threshold": "2.500000", 
            "u_substitution": "0", 
            "verbose_logs": "0"
        }, 
        "read_count": 7, 
        "reads_per_channel_dist": [
            {
                "channel": 33, 
                "count": 1
            }, 
            {
                "channel": 42, 
                "count": 1
            }, 
            {
                "channel": 55, 
                "count": 1
            }, 
            {
                "channel": 75, 
                "count": 1
            }, 
            {
                "channel": 95, 
                "count": 1
            }, 
            {
                "channel": 99, 
                "count": 1
            }, 
            {
                "channel": 109, 
                "count": 1
            }
        ], 
        "run_id": "6a94eb36563a9975b38e81c41534f6590c24867c", 
        "segment_duration": 60, 
        "segment_number": 3, 
        "segment_type": "guppy-acquisition", 
        "software": {
            "analysis": "1d_basecalling", 
            "name": "guppy-basecalling", 
            "version": "6.1.2+e0556ff"
        }, 
        "tracking_id": {
            "asic_id": "3506339664", 
            "asic_id_eeprom": "5924305", 
            "asic_temp": "35.086052", 
            "asic_version": "IA02D", 
            "auto_update": "0", 
            "auto_update_source": "https://mirror.oxfordnanoportal.com/software/MinKNOW/", 
            "bream_is_standard": "0", 
            "configuration_version": "4.2.11", 
            "device_id": "MN34875", 
            "device_type": "minion", 
            "distribution_status": "stable", 
            "distribution_version": "21.02.1", 
            "exp_script_name": "sequencing_MIN106_DNA:FLO-FLG001:SQK-LSK110", 
            "exp_script_purpose": "sequencing_run", 
            "exp_start_time": "2021-08-03T19:16:19Z", 
            "flongle_adapter_id": "FA-00993", 
            "flow_cell_id": "AHC239", 
            "flow_cell_product_code": "FLO-FLG001", 
            "guppy_version": "4.3.4+ecb2805", 
            "heatsink_temp": "35.003906", 
            "hostname": "mckennalab", 
            "installation_type": "nc", 
            "local_firmware_file": "1", 
            "msg_id": "f2d1c851-adaf-4c84-9ed6-07888ee2d34d", 
            "operating_system": "ubuntu 20.04", 
            "protocol_group_id": "FEX_120", 
            "protocol_run_id": "bdcd3859-b6c1-400c-b4e6-2bc287733f75", 
            "protocol_start_time": "", 
            "protocols_version": "6.1.10", 
            "run_id": "6a94eb36563a9975b38e81c41534f6590c24867c", 
            "sample_id": "96_v2", 
            "time_stamp": "2022-06-11T02:04:21Z", 
            "usb_config": "MinION_fx3_1.1.1_ONT#MinION_fpga_1.1.0#bulk#Auto", 
            "version": "4.2.5"
        }
    }, 
    {
        "aggregation": "segment", 
        "analysis_id": "e2692740-f853-4dfa-99c2-3309f56cdfdb", 
        "basecall_1d": {
            "exit_status_dist": {
                "pass": 6
            }, 
            "qscore_dist_temp": [
                {
                    "count": 2, 
                    "mean_qscore": 11.5
                }, 
                {
                    "count": 1, 
                    "mean_qscore": 12.5
                }, 
                {
                    "count": 2, 
                    "mean_qscore": 13.5
                }, 
                {
                    "count": 1, 
                    "mean_qscore": 14.0
                }
            ], 
            "qscore_sum_temp": {
                "count": 6, 
                "mean": 12.98987865447998, 
                "sum": 77.93927001953125
            }, 
            "read_len_events_sum_temp": 155440, 
            "seq_len_bases_dist_temp": [
                {
                    "count": 6, 
                    "length": 0.0
                }
            ], 
            "seq_len_bases_sum_temp": 6, 
            "seq_len_events_dist_temp": [
                {
                    "count": 1, 
                    "length": 17000.0
                }, 
                {
                    "count": 1, 
                    "length": 20000.0
                }, 
                {
                    "count": 2, 
                    "length": 25000.0
                }, 
                {
                    "count": 1, 
                    "length": 27000.0
                }, 
                {
                    "count": 1, 
                    "length": 38000.0
                }
            ], 
            "speed_bases_per_second_dist_temp": [
                {
                    "count": 6, 
                    "speed": 1.0
                }
            ], 
            "strand_median_pa": {
                "count": 6, 
                "mean": 70.73360443115234, 
                "sum": 424.4016418457031
            }, 
            "strand_sd_pa": {
                "count": 6, 
                "mean": 9.224841117858887, 
                "sum": 55.34904861450195
            }
        }, 
        "channel_count": 5, 
        "context_tags": {
            "barcoding_enabled": "0", 
            "experiment_duration_set": "1680", 
            "experiment_type": "genomic_dna", 
            "local_basecalling": "0", 
            "package": "bream4", 
            "package_version": "6.1.10", 
            "sample_frequency": "4000", 
            "sequencing_kit": "sqk-lsk110"
        }, 
        "latest_run_time": 14354.5966796875, 
        "levels_sums": {
            "count": 6, 
            "mean": 193.84765625, 
            "open_pore_level_sum": 1163.0859375
        }, 
        "opts": {
            "adapter_pt_range_scale": "5.200000", 
            "additional_lamp_context_bases": "2", 
            "align_ref": "", 
            "align_type": "auto", 
            "allow_inferior_barcodes": "0", 
            "as_cpu_threads_per_scaler": "2", 
            "as_gpu_runners_per_device": "2", 
            "as_model_file": "", 
            "as_num_scalers": "4", 
            "as_reads_per_runner": "32", 
            "bam_methylation_threshold": "5.000000", 
            "bam_out": "0", 
            "barcode_kits": "", 
            "barcode_nested_output_folder": "0", 
            "beam_cut": "100.000000", 
            "beam_width": "32", 
            "bed_file": "", 
            "builtin_scripts": "1", 
            "calib_detect": "0", 
            "calib_max_sequence_length": "3800", 
            "calib_min_coverage": "0.600000", 
            "calib_min_sequence_length": "3000", 
            "calib_reference": "lambda_3.6kb.fasta", 
            "chunk_size": "2000", 
            "chunks_per_caller": "10000", 
            "chunks_per_runner": "208", 
            "client_id": "-1", 
            "compress_fastq": "1", 
            "cpu_threads_per_caller": "4", 
            "detect_adapter": "0", 
            "detect_barcodes": "0", 
            "detect_mid_strand_adapter": "0", 
            "detect_mid_strand_barcodes": "0", 
            "detect_primer": "0", 
            "device": "cuda:0", 
            "disable_pings": "0", 
            "disable_qscore_filtering": "0", 
            "dmean_threshold": "1.000000", 
            "dmean_win_size": "2", 
            "do_read_splitting": "0", 
            "duplex_window_size_max": "1000", 
            "duplex_window_size_min": "200", 
            "end_gap1": "40", 
            "end_gap2": "40", 
            "extend_gap1": "40", 
            "extend_gap2": "160", 
            "fast5_out": "0", 
            "flowcell": "", 
            "front_window_size": "150", 
            "gpu_runners_per_device": "12", 
            "high_priority_threshold": "10", 
            "index": "0", 
            "input_file_list": "", 
            "int8_mode": "0", 
            "jump_threshold": "1.000000", 
            "kernel_path": "", 
            "kit": "", 
            "lamp_kit": "", 
            "log_speed_frequency": "0", 
            "max_queued_reads": "2000", 
            "max_read_split_depth": "2", 
            "max_search_len": "1000", 
            "medium_priority_threshold": "4", 
            "min_length_lamp_context": "30", 
            "min_length_lamp_target": "70", 
            "min_qscore": "10.000000", 
            "min_score_adapter": "60.000000", 
            "min_score_adapter_mid": "50.000000", 
            "min_score_barcode_front": "60.000000", 
            "min_score_barcode_mask": "40.000000", 
            "min_score_barcode_mid": "50.000000", 
            "min_score_barcode_rear": "60.000000", 
            "min_score_lamp": "80.000000", 
            "min_score_lamp_mask": "50.000000", 
            "min_score_lamp_target": "50.000000", 
            "min_score_primer": "60.000000", 
            "min_score_read_splitting": "70.000000", 
            "model_file": "template_r9.4.1_450bps_sup.jsn", 
            "moves_out": "0", 
            "nested_output_folder": "0", 
            "noisiest_section_scaling_max_size": "8000", 
            "num_alignment_threads": "4", 
            "num_barcode_threads": "4", 
            "num_barcoding_buffers": "24", 
            "num_base_mod_threads": "2", 
            "num_callers": "4", 
            "num_extra_bases_trim": "0", 
            "num_mid_barcoding_buffers": "96", 
            "num_read_splitting_buffers": "16", 
            "num_read_splitting_threads": "4", 
            "num_reads_per_barcoding_buffer": "4", 
            "open_gap1": "40", 
            "open_gap2": "160", 
            "overlap": "100", 
            "override_scaling": "0", 
            "ping_segment_duration": "60", 
            "ping_url": "https://ping.oxfordnanoportal.com/basecall", 
            "post_out": "0", 
            "print_workflows": "0", 
            "progress_stats_frequency": "-1.000000", 
            "pt_median_offset": "2.500000", 
            "pt_minimum_read_start_index": "30", 
            "pt_required_adapter_drop": "30.000000", 
            "pt_scaling": "0", 
            "qscore_offset": "0.349800", 
            "qscore_scale": "0.972200", 
            "quiet": "0", 
            "read_batch_size": "4000", 
            "read_id_list": "", 
            "read_splitting_arrangement_files": "", 
            "read_splitting_score_matrix_filename": "", 
            "rear_window_size": "150", 
            "records_per_fastq": "0", 
            "recursive": "0", 
            "remora_models": "", 
            "require_barcodes_both_ends": "0", 
            "resume": "0", 
            "reverse_sequence": "0", 
            "sample_sheet": "", 
            "scaling_mad": "1.000000", 
            "scaling_med": "0.000000", 
            "score_matrix_filename": "5x5_mismatch_matrix.txt", 
            "start_gap1": "40", 
            "start_gap2": "40", 
            "stay_penalty": "1.000000", 
            "temp_bias": "1.000000", 
            "temp_weight": "1.000000", 
            "trace_categories_logs": "", 
            "trace_domains_config": "", 
            "trim_adapters": "0", 
            "trim_barcodes": "0", 
            "trim_min_events": "3", 
            "trim_primers": "0", 
            "trim_strategy": "dna", 
            "trim_threshold": "2.500000", 
            "u_substitution": "0", 
            "verbose_logs": "0"
        }, 
        "read_count": 6, 
        "reads_per_channel_dist": [
            {
                "channel": 54, 
                "count": 1
            }, 
            {
                "channel": 56, 
                "count": 1
            }, 
            {
                "channel": 62, 
                "count": 1
            }, 
            {
                "channel": 66, 
                "count": 1
            }, 
            {
                "channel": 111, 
                "count": 2
            }
        ], 
        "run_id": "6a94eb36563a9975b38e81c41534f6590c24867c", 
        "segment_duration": 60, 
        "segment_number": 4, 
        "segment_type": "guppy-acquisition", 
        "software": {
            "analysis": "1d_basecalling", 
            "name": "guppy-basecalling", 
            "version": "6.1.2+e0556ff"
        }, 
        "tracking_id": {
            "asic_id": "3506339664", 
            "asic_id_eeprom": "5924305", 
            "asic_temp": "35.086052", 
            "asic_version": "IA02D", 
            "auto_update": "0", 
            "auto_update_source": "https://mirror.oxfordnanoportal.com/software/MinKNOW/", 
            "bream_is_standard": "0", 
            "configuration_version": "4.2.11", 
            "device_id": "MN34875", 
            "device_type": "minion", 
            "distribution_status": "stable", 
            "distribution_version": "21.02.1", 
            "exp_script_name": "sequencing_MIN106_DNA:FLO-FLG001:SQK-LSK110", 
            "exp_script_purpose": "sequencing_run", 
            "exp_start_time": "2021-08-03T19:16:19Z", 
            "flongle_adapter_id": "FA-00993", 
            "flow_cell_id": "AHC239", 
            "flow_cell_product_code": "FLO-FLG001", 
            "guppy_version": "4.3.4+ecb2805", 
            "heatsink_temp": "35.003906", 
            "hostname": "mckennalab", 
            "installation_type": "nc", 
            "local_firmware_file": "1", 
            "msg_id": "0c5f339f-21f4-40fb-be4d-972973e6e6b8", 
            "operating_system": "ubuntu 20.04", 
            "protocol_group_id": "FEX_120", 
            "protocol_run_id": "bdcd3859-b6c1-400c-b4e6-2bc287733f75", 
            "protocol_start_time": "", 
            "protocols_version": "6.1.10", 
            "run_id": "6a94eb36563a9975b38e81c41534f6590c24867c", 
            "sample_id": "96_v2", 
            "time_stamp": "2022-06-11T02:04:21Z", 
            "usb_config": "MinION_fx3_1.1.1_ONT#MinION_fpga_1.1.0#bulk#Auto", 
            "version": "4.2.5"
        }
    }, 
    {
        "aggregation": "segment", 
        "analysis_id": "e2692740-f853-4dfa-99c2-3309f56cdfdb", 
        "basecall_1d": {
            "exit_status_dist": {
                "pass": 1
            }, 
            "qscore_dist_temp": [
                {
                    "count": 1, 
                    "mean_qscore": 12.0
                }
            ], 
            "qscore_sum_temp": {
                "count": 1, 
                "mean": 12.462763786315918, 
                "sum": 12.462763786315918
            }, 
            "read_len_events_sum_temp": 17169, 
            "seq_len_bases_dist_temp": [
                {
                    "count": 1, 
                    "length": 0.0
                }
            ], 
            "seq_len_bases_sum_temp": 1, 
            "seq_len_events_dist_temp": [
                {
                    "count": 1, 
                    "length": 17000.0
                }
            ], 
            "speed_bases_per_second_dist_temp": [
                {
                    "count": 1, 
                    "speed": 1.0
                }
            ], 
            "strand_median_pa": {
                "count": 1, 
                "mean": 73.0324478149414, 
                "sum": 73.0324478149414
            }, 
            "strand_sd_pa": {
                "count": 1, 
                "mean": 9.195368766784668, 
                "sum": 9.195368766784668
            }
        }, 
        "channel_count": 1, 
        "context_tags": {
            "barcoding_enabled": "0", 
            "experiment_duration_set": "1680", 
            "experiment_type": "genomic_dna", 
            "local_basecalling": "0", 
            "package": "bream4", 
            "package_version": "6.1.10", 
            "sample_frequency": "4000", 
            "sequencing_kit": "sqk-lsk110"
        }, 
        "latest_run_time": 14985.69921875, 
        "levels_sums": {
            "count": 1, 
            "mean": 195.72210693359375, 
            "open_pore_level_sum": 195.72210693359375
        }, 
        "opts": {
            "adapter_pt_range_scale": "5.200000", 
            "additional_lamp_context_bases": "2", 
            "align_ref": "", 
            "align_type": "auto", 
            "allow_inferior_barcodes": "0", 
            "as_cpu_threads_per_scaler": "2", 
            "as_gpu_runners_per_device": "2", 
            "as_model_file": "", 
            "as_num_scalers": "4", 
            "as_reads_per_runner": "32", 
            "bam_methylation_threshold": "5.000000", 
            "bam_out": "0", 
            "barcode_kits": "", 
            "barcode_nested_output_folder": "0", 
            "beam_cut": "100.000000", 
            "beam_width": "32", 
            "bed_file": "", 
            "builtin_scripts": "1", 
            "calib_detect": "0", 
            "calib_max_sequence_length": "3800", 
            "calib_min_coverage": "0.600000", 
            "calib_min_sequence_length": "3000", 
            "calib_reference": "lambda_3.6kb.fasta", 
            "chunk_size": "2000", 
            "chunks_per_caller": "10000", 
            "chunks_per_runner": "208", 
            "client_id": "-1", 
            "compress_fastq": "1", 
            "cpu_threads_per_caller": "4", 
            "detect_adapter": "0", 
            "detect_barcodes": "0", 
            "detect_mid_strand_adapter": "0", 
            "detect_mid_strand_barcodes": "0", 
            "detect_primer": "0", 
            "device": "cuda:0", 
            "disable_pings": "0", 
            "disable_qscore_filtering": "0", 
            "dmean_threshold": "1.000000", 
            "dmean_win_size": "2", 
            "do_read_splitting": "0", 
            "duplex_window_size_max": "1000", 
            "duplex_window_size_min": "200", 
            "end_gap1": "40", 
            "end_gap2": "40", 
            "extend_gap1": "40", 
            "extend_gap2": "160", 
            "fast5_out": "0", 
            "flowcell": "", 
            "front_window_size": "150", 
            "gpu_runners_per_device": "12", 
            "high_priority_threshold": "10", 
            "index": "0", 
            "input_file_list": "", 
            "int8_mode": "0", 
            "jump_threshold": "1.000000", 
            "kernel_path": "", 
            "kit": "", 
            "lamp_kit": "", 
            "log_speed_frequency": "0", 
            "max_queued_reads": "2000", 
            "max_read_split_depth": "2", 
            "max_search_len": "1000", 
            "medium_priority_threshold": "4", 
            "min_length_lamp_context": "30", 
            "min_length_lamp_target": "70", 
            "min_qscore": "10.000000", 
            "min_score_adapter": "60.000000", 
            "min_score_adapter_mid": "50.000000", 
            "min_score_barcode_front": "60.000000", 
            "min_score_barcode_mask": "40.000000", 
            "min_score_barcode_mid": "50.000000", 
            "min_score_barcode_rear": "60.000000", 
            "min_score_lamp": "80.000000", 
            "min_score_lamp_mask": "50.000000", 
            "min_score_lamp_target": "50.000000", 
            "min_score_primer": "60.000000", 
            "min_score_read_splitting": "70.000000", 
            "model_file": "template_r9.4.1_450bps_sup.jsn", 
            "moves_out": "0", 
            "nested_output_folder": "0", 
            "noisiest_section_scaling_max_size": "8000", 
            "num_alignment_threads": "4", 
            "num_barcode_threads": "4", 
            "num_barcoding_buffers": "24", 
            "num_base_mod_threads": "2", 
            "num_callers": "4", 
            "num_extra_bases_trim": "0", 
            "num_mid_barcoding_buffers": "96", 
            "num_read_splitting_buffers": "16", 
            "num_read_splitting_threads": "4", 
            "num_reads_per_barcoding_buffer": "4", 
            "open_gap1": "40", 
            "open_gap2": "160", 
            "overlap": "100", 
            "override_scaling": "0", 
            "ping_segment_duration": "60", 
            "ping_url": "https://ping.oxfordnanoportal.com/basecall", 
            "post_out": "0", 
            "print_workflows": "0", 
            "progress_stats_frequency": "-1.000000", 
            "pt_median_offset": "2.500000", 
            "pt_minimum_read_start_index": "30", 
            "pt_required_adapter_drop": "30.000000", 
            "pt_scaling": "0", 
            "qscore_offset": "0.349800", 
            "qscore_scale": "0.972200", 
            "quiet": "0", 
            "read_batch_size": "4000", 
            "read_id_list": "", 
            "read_splitting_arrangement_files": "", 
            "read_splitting_score_matrix_filename": "", 
            "rear_window_size": "150", 
            "records_per_fastq": "0", 
            "recursive": "0", 
            "remora_models": "", 
            "require_barcodes_both_ends": "0", 
            "resume": "0", 
            "reverse_sequence": "0", 
            "sample_sheet": "", 
            "scaling_mad": "1.000000", 
            "scaling_med": "0.000000", 
            "score_matrix_filename": "5x5_mismatch_matrix.txt", 
            "start_gap1": "40", 
            "start_gap2": "40", 
            "stay_penalty": "1.000000", 
            "temp_bias": "1.000000", 
            "temp_weight": "1.000000", 
            "trace_categories_logs": "", 
            "trace_domains_config": "", 
            "trim_adapters": "0", 
            "trim_barcodes": "0", 
            "trim_min_events": "3", 
            "trim_primers": "0", 
            "trim_strategy": "dna", 
            "trim_threshold": "2.500000", 
            "u_substitution": "0", 
            "verbose_logs": "0"
        }, 
        "read_count": 1, 
        "reads_per_channel_dist": [
            {
                "channel": 17, 
                "count": 1
            }
        ], 
        "run_id": "6a94eb36563a9975b38e81c41534f6590c24867c", 
        "segment_duration": 60, 
        "segment_number": 5, 
        "segment_type": "guppy-acquisition", 
        "software": {
            "analysis": "1d_basecalling", 
            "name": "guppy-basecalling", 
            "version": "6.1.2+e0556ff"
        }, 
        "tracking_id": {
            "asic_id": "3506339664", 
            "asic_id_eeprom": "5924305", 
            "asic_temp": "35.086052", 
            "asic_version": "IA02D", 
            "auto_update": "0", 
            "auto_update_source": "https://mirror.oxfordnanoportal.com/software/MinKNOW/", 
            "bream_is_standard": "0", 
            "configuration_version": "4.2.11", 
            "device_id": "MN34875", 
            "device_type": "minion", 
            "distribution_status": "stable", 
            "distribution_version": "21.02.1", 
            "exp_script_name": "sequencing_MIN106_DNA:FLO-FLG001:SQK-LSK110", 
            "exp_script_purpose": "sequencing_run", 
            "exp_start_time": "2021-08-03T19:16:19Z", 
            "flongle_adapter_id": "FA-00993", 
            "flow_cell_id": "AHC239", 
            "flow_cell_product_code": "FLO-FLG001", 
            "guppy_version": "4.3.4+ecb2805", 
            "heatsink_temp": "35.003906", 
            "hostname": "mckennalab", 
            "installation_type": "nc", 
            "local_firmware_file": "1", 
            "msg_id": "7e2b0ffd-aa3d-4bd3-bd65-4dee16379f0c", 
            "operating_system": "ubuntu 20.04", 
            "protocol_group_id": "FEX_120", 
            "protocol_run_id": "bdcd3859-b6c1-400c-b4e6-2bc287733f75", 
            "protocol_start_time": "", 
            "protocols_version": "6.1.10", 
            "run_id": "6a94eb36563a9975b38e81c41534f6590c24867c", 
            "sample_id": "96_v2", 
            "time_stamp": "2022-06-11T02:04:21Z", 
            "usb_config": "MinION_fx3_1.1.1_ONT#MinION_fpga_1.1.0#bulk#Auto", 
            "version": "4.2.5"
        }
    }, 
    {
        "aggregation": "segment", 
        "analysis_id": "e2692740-f853-4dfa-99c2-3309f56cdfdb", 
        "basecall_1d": {
            "exit_status_dist": {
                "pass": 5
            }, 
            "qscore_dist_temp": [
                {
                    "count": 1, 
                    "mean_qscore": 11.5
                }, 
                {
                    "count": 2, 
                    "mean_qscore": 12.5
                }, 
                {
                    "count": 1, 
                    "mean_qscore": 13.0
                }, 
                {
                    "count": 1, 
                    "mean_qscore": 13.5
                }
            ], 
            "qscore_sum_temp": {
                "count": 5, 
                "mean": 12.81437873840332, 
                "sum": 64.07189178466797
            }, 
            "read_len_events_sum_temp": 119132, 
            "seq_len_bases_dist_temp": [
                {
                    "count": 5, 
                    "length": 0.0
                }
            ], 
            "seq_len_bases_sum_temp": 5, 
            "seq_len_events_dist_temp": [
                {
                    "count": 1, 
                    "length": 16000.0
                }, 
                {
                    "count": 1, 
                    "length": 17000.0
                }, 
                {
                    "count": 2, 
                    "length": 18000.0
                }, 
                {
                    "count": 1, 
                    "length": 47000.0
                }
            ], 
            "speed_bases_per_second_dist_temp": [
                {
                    "count": 5, 
                    "speed": 1.0
                }
            ], 
            "strand_median_pa": {
                "count": 5, 
                "mean": 69.81407165527344, 
                "sum": 349.0703430175781
            }, 
            "strand_sd_pa": {
                "count": 5, 
                "mean": 9.195368766784668, 
                "sum": 45.976844787597656
            }
        }, 
        "channel_count": 5, 
        "context_tags": {
            "barcoding_enabled": "0", 
            "experiment_duration_set": "1680", 
            "experiment_type": "genomic_dna", 
            "local_basecalling": "0", 
            "package": "bream4", 
            "package_version": "6.1.10", 
            "sample_frequency": "4000", 
            "sequencing_kit": "sqk-lsk110"
        }, 
        "latest_run_time": 20825.4375, 
        "levels_sums": {
            "count": 5, 
            "mean": 193.79000854492188, 
            "open_pore_level_sum": 968.9500732421875
        }, 
        "opts": {
            "adapter_pt_range_scale": "5.200000", 
            "additional_lamp_context_bases": "2", 
            "align_ref": "", 
            "align_type": "auto", 
            "allow_inferior_barcodes": "0", 
            "as_cpu_threads_per_scaler": "2", 
            "as_gpu_runners_per_device": "2", 
            "as_model_file": "", 
            "as_num_scalers": "4", 
            "as_reads_per_runner": "32", 
            "bam_methylation_threshold": "5.000000", 
            "bam_out": "0", 
            "barcode_kits": "", 
            "barcode_nested_output_folder": "0", 
            "beam_cut": "100.000000", 
            "beam_width": "32", 
            "bed_file": "", 
            "builtin_scripts": "1", 
            "calib_detect": "0", 
            "calib_max_sequence_length": "3800", 
            "calib_min_coverage": "0.600000", 
            "calib_min_sequence_length": "3000", 
            "calib_reference": "lambda_3.6kb.fasta", 
            "chunk_size": "2000", 
            "chunks_per_caller": "10000", 
            "chunks_per_runner": "208", 
            "client_id": "-1", 
            "compress_fastq": "1", 
            "cpu_threads_per_caller": "4", 
            "detect_adapter": "0", 
            "detect_barcodes": "0", 
            "detect_mid_strand_adapter": "0", 
            "detect_mid_strand_barcodes": "0", 
            "detect_primer": "0", 
            "device": "cuda:0", 
            "disable_pings": "0", 
            "disable_qscore_filtering": "0", 
            "dmean_threshold": "1.000000", 
            "dmean_win_size": "2", 
            "do_read_splitting": "0", 
            "duplex_window_size_max": "1000", 
            "duplex_window_size_min": "200", 
            "end_gap1": "40", 
            "end_gap2": "40", 
            "extend_gap1": "40", 
            "extend_gap2": "160", 
            "fast5_out": "0", 
            "flowcell": "", 
            "front_window_size": "150", 
            "gpu_runners_per_device": "12", 
            "high_priority_threshold": "10", 
            "index": "0", 
            "input_file_list": "", 
            "int8_mode": "0", 
            "jump_threshold": "1.000000", 
            "kernel_path": "", 
            "kit": "", 
            "lamp_kit": "", 
            "log_speed_frequency": "0", 
            "max_queued_reads": "2000", 
            "max_read_split_depth": "2", 
            "max_search_len": "1000", 
            "medium_priority_threshold": "4", 
            "min_length_lamp_context": "30", 
            "min_length_lamp_target": "70", 
            "min_qscore": "10.000000", 
            "min_score_adapter": "60.000000", 
            "min_score_adapter_mid": "50.000000", 
            "min_score_barcode_front": "60.000000", 
            "min_score_barcode_mask": "40.000000", 
            "min_score_barcode_mid": "50.000000", 
            "min_score_barcode_rear": "60.000000", 
            "min_score_lamp": "80.000000", 
            "min_score_lamp_mask": "50.000000", 
            "min_score_lamp_target": "50.000000", 
            "min_score_primer": "60.000000", 
            "min_score_read_splitting": "70.000000", 
            "model_file": "template_r9.4.1_450bps_sup.jsn", 
            "moves_out": "0", 
            "nested_output_folder": "0", 
            "noisiest_section_scaling_max_size": "8000", 
            "num_alignment_threads": "4", 
            "num_barcode_threads": "4", 
            "num_barcoding_buffers": "24", 
            "num_base_mod_threads": "2", 
            "num_callers": "4", 
            "num_extra_bases_trim": "0", 
            "num_mid_barcoding_buffers": "96", 
            "num_read_splitting_buffers": "16", 
            "num_read_splitting_threads": "4", 
            "num_reads_per_barcoding_buffer": "4", 
            "open_gap1": "40", 
            "open_gap2": "160", 
            "overlap": "100", 
            "override_scaling": "0", 
            "ping_segment_duration": "60", 
            "ping_url": "https://ping.oxfordnanoportal.com/basecall", 
            "post_out": "0", 
            "print_workflows": "0", 
            "progress_stats_frequency": "-1.000000", 
            "pt_median_offset": "2.500000", 
            "pt_minimum_read_start_index": "30", 
            "pt_required_adapter_drop": "30.000000", 
            "pt_scaling": "0", 
            "qscore_offset": "0.349800", 
            "qscore_scale": "0.972200", 
            "quiet": "0", 
            "read_batch_size": "4000", 
            "read_id_list": "", 
            "read_splitting_arrangement_files": "", 
            "read_splitting_score_matrix_filename": "", 
            "rear_window_size": "150", 
            "records_per_fastq": "0", 
            "recursive": "0", 
            "remora_models": "", 
            "require_barcodes_both_ends": "0", 
            "resume": "0", 
            "reverse_sequence": "0", 
            "sample_sheet": "", 
            "scaling_mad": "1.000000", 
            "scaling_med": "0.000000", 
            "score_matrix_filename": "5x5_mismatch_matrix.txt", 
            "start_gap1": "40", 
            "start_gap2": "40", 
            "stay_penalty": "1.000000", 
            "temp_bias": "1.000000", 
            "temp_weight": "1.000000", 
            "trace_categories_logs": "", 
            "trace_domains_config": "", 
            "trim_adapters": "0", 
            "trim_barcodes": "0", 
            "trim_min_events": "3", 
            "trim_primers": "0", 
            "trim_strategy": "dna", 
            "trim_threshold": "2.500000", 
            "u_substitution": "0", 
            "verbose_logs": "0"
        }, 
        "read_count": 5, 
        "reads_per_channel_dist": [
            {
                "channel": 23, 
                "count": 1
            }, 
            {
                "channel": 27, 
                "count": 1
            }, 
            {
                "channel": 28, 
                "count": 1
            }, 
            {
                "channel": 49, 
                "count": 1
            }, 
            {
                "channel": 112, 
                "count": 1
            }
        ], 
        "run_id": "6a94eb36563a9975b38e81c41534f6590c24867c", 
        "segment_duration": 60, 
        "segment_number": 6, 
        "segment_type": "guppy-acquisition", 
        "software": {
            "analysis": "1d_basecalling", 
            "name": "guppy-basecalling", 
            "version": "6.1.2+e0556ff"
        }, 
        "tracking_id": {
            "asic_id": "3506339664", 
            "asic_id_eeprom": "5924305", 
            "asic_temp": "35.086052", 
            "asic_version": "IA02D", 
            "auto_update": "0", 
            "auto_update_source": "https://mirror.oxfordnanoportal.com/software/MinKNOW/", 
            "bream_is_standard": "0", 
            "configuration_version": "4.2.11", 
            "device_id": "MN34875", 
            "device_type": "minion", 
            "distribution_status": "stable", 
            "distribution_version": "21.02.1", 
            "exp_script_name": "sequencing_MIN106_DNA:FLO-FLG001:SQK-LSK110", 
            "exp_script_purpose": "sequencing_run", 
            "exp_start_time": "2021-08-03T19:16:19Z", 
            "flongle_adapter_id": "FA-00993", 
            "flow_cell_id": "AHC239", 
            "flow_cell_product_code": "FLO-FLG001", 
            "guppy_version": "4.3.4+ecb2805", 
            "heatsink_temp": "35.003906", 
            "hostname": "mckennalab", 
            "installation_type": "nc", 
            "local_firmware_file": "1", 
            "msg_id": "9e0f6385-541c-4ab8-9e77-8356b41dd66d", 
            "operating_system": "ubuntu 20.04", 
            "protocol_group_id": "FEX_120", 
            "protocol_run_id": "bdcd3859-b6c1-400c-b4e6-2bc287733f75", 
            "protocol_start_time": "", 
            "protocols_version": "6.1.10", 
            "run_id": "6a94eb36563a9975b38e81c41534f6590c24867c", 
            "sample_id": "96_v2", 
            "time_stamp": "2022-06-11T02:04:21Z", 
            "usb_config": "MinION_fx3_1.1.1_ONT#MinION_fpga_1.1.0#bulk#Auto", 
            "version": "4.2.5"
        }
    }, 
    {
        "aggregation": "segment", 
        "analysis_id": "e2692740-f853-4dfa-99c2-3309f56cdfdb", 
        "basecall_1d": {
            "exit_status_dist": {
                "pass": 1
            }, 
            "qscore_dist_temp": [
                {
                    "count": 1, 
                    "mean_qscore": 10.5
                }
            ], 
            "qscore_sum_temp": {
                "count": 1, 
                "mean": 10.597357749938965, 
                "sum": 10.597357749938965
            }, 
            "read_len_events_sum_temp": 21182, 
            "seq_len_bases_dist_temp": [
                {
                    "count": 1, 
                    "length": 0.0
                }
            ], 
            "seq_len_bases_sum_temp": 1, 
            "seq_len_events_dist_temp": [
                {
                    "count": 1, 
                    "length": 21000.0
                }
            ], 
            "speed_bases_per_second_dist_temp": [
                {
                    "count": 1, 
                    "speed": 1.0
                }
            ], 
            "strand_median_pa": {
                "count": 1, 
                "mean": 77.63013458251953, 
                "sum": 77.63013458251953
            }, 
            "strand_sd_pa": {
                "count": 1, 
                "mean": 9.37220287322998, 
                "sum": 9.37220287322998
            }
        }, 
        "channel_count": 1, 
        "context_tags": {
            "barcoding_enabled": "0", 
            "experiment_duration_set": "1680", 
            "experiment_type": "genomic_dna", 
            "local_basecalling": "0", 
            "package": "bream4", 
            "package_version": "6.1.10", 
            "sample_frequency": "4000", 
            "sequencing_kit": "sqk-lsk110"
        }, 
        "latest_run_time": 23201.220703125, 
        "levels_sums": {
            "count": 1, 
            "mean": 209.0736541748047, 
            "open_pore_level_sum": 209.0736541748047
        }, 
        "opts": {
            "adapter_pt_range_scale": "5.200000", 
            "additional_lamp_context_bases": "2", 
            "align_ref": "", 
            "align_type": "auto", 
            "allow_inferior_barcodes": "0", 
            "as_cpu_threads_per_scaler": "2", 
            "as_gpu_runners_per_device": "2", 
            "as_model_file": "", 
            "as_num_scalers": "4", 
            "as_reads_per_runner": "32", 
            "bam_methylation_threshold": "5.000000", 
            "bam_out": "0", 
            "barcode_kits": "", 
            "barcode_nested_output_folder": "0", 
            "beam_cut": "100.000000", 
            "beam_width": "32", 
            "bed_file": "", 
            "builtin_scripts": "1", 
            "calib_detect": "0", 
            "calib_max_sequence_length": "3800", 
            "calib_min_coverage": "0.600000", 
            "calib_min_sequence_length": "3000", 
            "calib_reference": "lambda_3.6kb.fasta", 
            "chunk_size": "2000", 
            "chunks_per_caller": "10000", 
            "chunks_per_runner": "208", 
            "client_id": "-1", 
            "compress_fastq": "1", 
            "cpu_threads_per_caller": "4", 
            "detect_adapter": "0", 
            "detect_barcodes": "0", 
            "detect_mid_strand_adapter": "0", 
            "detect_mid_strand_barcodes": "0", 
            "detect_primer": "0", 
            "device": "cuda:0", 
            "disable_pings": "0", 
            "disable_qscore_filtering": "0", 
            "dmean_threshold": "1.000000", 
            "dmean_win_size": "2", 
            "do_read_splitting": "0", 
            "duplex_window_size_max": "1000", 
            "duplex_window_size_min": "200", 
            "end_gap1": "40", 
            "end_gap2": "40", 
            "extend_gap1": "40", 
            "extend_gap2": "160", 
            "fast5_out": "0", 
            "flowcell": "", 
            "front_window_size": "150", 
            "gpu_runners_per_device": "12", 
            "high_priority_threshold": "10", 
            "index": "0", 
            "input_file_list": "", 
            "int8_mode": "0", 
            "jump_threshold": "1.000000", 
            "kernel_path": "", 
            "kit": "", 
            "lamp_kit": "", 
            "log_speed_frequency": "0", 
            "max_queued_reads": "2000", 
            "max_read_split_depth": "2", 
            "max_search_len": "1000", 
            "medium_priority_threshold": "4", 
            "min_length_lamp_context": "30", 
            "min_length_lamp_target": "70", 
            "min_qscore": "10.000000", 
            "min_score_adapter": "60.000000", 
            "min_score_adapter_mid": "50.000000", 
            "min_score_barcode_front": "60.000000", 
            "min_score_barcode_mask": "40.000000", 
            "min_score_barcode_mid": "50.000000", 
            "min_score_barcode_rear": "60.000000", 
            "min_score_lamp": "80.000000", 
            "min_score_lamp_mask": "50.000000", 
            "min_score_lamp_target": "50.000000", 
            "min_score_primer": "60.000000", 
            "min_score_read_splitting": "70.000000", 
            "model_file": "template_r9.4.1_450bps_sup.jsn", 
            "moves_out": "0", 
            "nested_output_folder": "0", 
            "noisiest_section_scaling_max_size": "8000", 
            "num_alignment_threads": "4", 
            "num_barcode_threads": "4", 
            "num_barcoding_buffers": "24", 
            "num_base_mod_threads": "2", 
            "num_callers": "4", 
            "num_extra_bases_trim": "0", 
            "num_mid_barcoding_buffers": "96", 
            "num_read_splitting_buffers": "16", 
            "num_read_splitting_threads": "4", 
            "num_reads_per_barcoding_buffer": "4", 
            "open_gap1": "40", 
            "open_gap2": "160", 
            "overlap": "100", 
            "override_scaling": "0", 
            "ping_segment_duration": "60", 
            "ping_url": "https://ping.oxfordnanoportal.com/basecall", 
            "post_out": "0", 
            "print_workflows": "0", 
            "progress_stats_frequency": "-1.000000", 
            "pt_median_offset": "2.500000", 
            "pt_minimum_read_start_index": "30", 
            "pt_required_adapter_drop": "30.000000", 
            "pt_scaling": "0", 
            "qscore_offset": "0.349800", 
            "qscore_scale": "0.972200", 
            "quiet": "0", 
            "read_batch_size": "4000", 
            "read_id_list": "", 
            "read_splitting_arrangement_files": "", 
            "read_splitting_score_matrix_filename": "", 
            "rear_window_size": "150", 
            "records_per_fastq": "0", 
            "recursive": "0", 
            "remora_models": "", 
            "require_barcodes_both_ends": "0", 
            "resume": "0", 
            "reverse_sequence": "0", 
            "sample_sheet": "", 
            "scaling_mad": "1.000000", 
            "scaling_med": "0.000000", 
            "score_matrix_filename": "5x5_mismatch_matrix.txt", 
            "start_gap1": "40", 
            "start_gap2": "40", 
            "stay_penalty": "1.000000", 
            "temp_bias": "1.000000", 
            "temp_weight": "1.000000", 
            "trace_categories_logs": "", 
            "trace_domains_config": "", 
            "trim_adapters": "0", 
            "trim_barcodes": "0", 
            "trim_min_events": "3", 
            "trim_primers": "0", 
            "trim_strategy": "dna", 
            "trim_threshold": "2.500000", 
            "u_substitution": "0", 
            "verbose_logs": "0"
        }, 
        "read_count": 1, 
        "reads_per_channel_dist": [
            {
                "channel": 26, 
                "count": 1
            }
        ], 
        "run_id": "6a94eb36563a9975b38e81c41534f6590c24867c", 
        "segment_duration": 60, 
        "segment_number": 7, 
        "segment_type": "guppy-acquisition", 
        "software": {
            "analysis": "1d_basecalling", 
            "name": "guppy-basecalling", 
            "version": "6.1.2+e0556ff"
        }, 
        "tracking_id": {
            "asic_id": "3506339664", 
            "asic_id_eeprom": "5924305", 
            "asic_temp": "35.086052", 
            "asic_version": "IA02D", 
            "auto_update": "0", 
            "auto_update_source": "https://mirror.oxfordnanoportal.com/software/MinKNOW/", 
            "bream_is_standard": "0", 
            "configuration_version": "4.2.11", 
            "device_id": "MN34875", 
            "device_type": "minion", 
            "distribution_status": "stable", 
            "distribution_version": "21.02.1", 
            "exp_script_name": "sequencing_MIN106_DNA:FLO-FLG001:SQK-LSK110", 
            "exp_script_purpose": "sequencing_run", 
            "exp_start_time": "2021-08-03T19:16:19Z", 
            "flongle_adapter_id": "FA-00993", 
            "flow_cell_id": "AHC239", 
            "flow_cell_product_code": "FLO-FLG001", 
            "guppy_version": "4.3.4+ecb2805", 
            "heatsink_temp": "35.003906", 
            "hostname": "mckennalab", 
            "installation_type": "nc", 
            "local_firmware_file": "1", 
            "msg_id": "edd7b8da-dda4-4f11-98c4-fdc3b4eb0d23", 
            "operating_system": "ubuntu 20.04", 
            "protocol_group_id": "FEX_120", 
            "protocol_run_id": "bdcd3859-b6c1-400c-b4e6-2bc287733f75", 
            "protocol_start_time": "", 
            "protocols_version": "6.1.10", 
            "run_id": "6a94eb36563a9975b38e81c41534f6590c24867c", 
            "sample_id": "96_v2", 
            "time_stamp": "2022-06-11T02:04:21Z", 
            "usb_config": "MinION_fx3_1.1.1_ONT#MinION_fpga_1.1.0#bulk#Auto", 
            "version": "4.2.5"
        }
    }, 
    {
        "aggregation": "segment", 
        "analysis_id": "e2692740-f853-4dfa-99c2-3309f56cdfdb", 
        "basecall_1d": {
            "exit_status_dist": {
                "pass": 5
            }, 
            "qscore_dist_temp": [
                {
                    "count": 1, 
                    "mean_qscore": 10.0
                }, 
                {
                    "count": 2, 
                    "mean_qscore": 11.5
                }, 
                {
                    "count": 1, 
                    "mean_qscore": 12.0
                }, 
                {
                    "count": 1, 
                    "mean_qscore": 13.5
                }
            ], 
            "qscore_sum_temp": {
                "count": 5, 
                "mean": 12.115621566772461, 
                "sum": 60.57810592651367
            }, 
            "read_len_events_sum_temp": 124991, 
            "seq_len_bases_dist_temp": [
                {
                    "count": 5, 
                    "length": 0.0
                }
            ], 
            "seq_len_bases_sum_temp": 5, 
            "seq_len_events_dist_temp": [
                {
                    "count": 1, 
                    "length": 22000.0
                }, 
                {
                    "count": 2, 
                    "length": 24000.0
                }, 
                {
                    "count": 1, 
                    "length": 25000.0
                }, 
                {
                    "count": 1, 
                    "length": 27000.0
                }
            ], 
            "speed_bases_per_second_dist_temp": [
                {
                    "count": 5, 
                    "speed": 1.0
                }
            ], 
            "strand_median_pa": {
                "count": 5, 
                "mean": 73.77515411376953, 
                "sum": 368.8757629394531
            }, 
            "strand_sd_pa": {
                "count": 5, 
                "mean": 9.301469802856445, 
                "sum": 46.507347106933594
            }
        }, 
        "channel_count": 5, 
        "context_tags": {
            "barcoding_enabled": "0", 
            "experiment_duration_set": "1680", 
            "experiment_type": "genomic_dna", 
            "local_basecalling": "0", 
            "package": "bream4", 
            "package_version": "6.1.10", 
            "sample_frequency": "4000", 
            "sequencing_kit": "sqk-lsk110"
        }, 
        "latest_run_time": 28299.900390625, 
        "levels_sums": {
            "count": 5, 
            "mean": 174.89393615722656, 
            "open_pore_level_sum": 874.46966552734375
        }, 
        "opts": {
            "adapter_pt_range_scale": "5.200000", 
            "additional_lamp_context_bases": "2", 
            "align_ref": "", 
            "align_type": "auto", 
            "allow_inferior_barcodes": "0", 
            "as_cpu_threads_per_scaler": "2", 
            "as_gpu_runners_per_device": "2", 
            "as_model_file": "", 
            "as_num_scalers": "4", 
            "as_reads_per_runner": "32", 
            "bam_methylation_threshold": "5.000000", 
            "bam_out": "0", 
            "barcode_kits": "", 
            "barcode_nested_output_folder": "0", 
            "beam_cut": "100.000000", 
            "beam_width": "32", 
            "bed_file": "", 
            "builtin_scripts": "1", 
            "calib_detect": "0", 
            "calib_max_sequence_length": "3800", 
            "calib_min_coverage": "0.600000", 
            "calib_min_sequence_length": "3000", 
            "calib_reference": "lambda_3.6kb.fasta", 
            "chunk_size": "2000", 
            "chunks_per_caller": "10000", 
            "chunks_per_runner": "208", 
            "client_id": "-1", 
            "compress_fastq": "1", 
            "cpu_threads_per_caller": "4", 
            "detect_adapter": "0", 
            "detect_barcodes": "0", 
            "detect_mid_strand_adapter": "0", 
            "detect_mid_strand_barcodes": "0", 
            "detect_primer": "0", 
            "device": "cuda:0", 
            "disable_pings": "0", 
            "disable_qscore_filtering": "0", 
            "dmean_threshold": "1.000000", 
            "dmean_win_size": "2", 
            "do_read_splitting": "0", 
            "duplex_window_size_max": "1000", 
            "duplex_window_size_min": "200", 
            "end_gap1": "40", 
            "end_gap2": "40", 
            "extend_gap1": "40", 
            "extend_gap2": "160", 
            "fast5_out": "0", 
            "flowcell": "", 
            "front_window_size": "150", 
            "gpu_runners_per_device": "12", 
            "high_priority_threshold": "10", 
            "index": "0", 
            "input_file_list": "", 
            "int8_mode": "0", 
            "jump_threshold": "1.000000", 
            "kernel_path": "", 
            "kit": "", 
            "lamp_kit": "", 
            "log_speed_frequency": "0", 
            "max_queued_reads": "2000", 
            "max_read_split_depth": "2", 
            "max_search_len": "1000", 
            "medium_priority_threshold": "4", 
            "min_length_lamp_context": "30", 
            "min_length_lamp_target": "70", 
            "min_qscore": "10.000000", 
            "min_score_adapter": "60.000000", 
            "min_score_adapter_mid": "50.000000", 
            "min_score_barcode_front": "60.000000", 
            "min_score_barcode_mask": "40.000000", 
            "min_score_barcode_mid": "50.000000", 
            "min_score_barcode_rear": "60.000000", 
            "min_score_lamp": "80.000000", 
            "min_score_lamp_mask": "50.000000", 
            "min_score_lamp_target": "50.000000", 
            "min_score_primer": "60.000000", 
            "min_score_read_splitting": "70.000000", 
            "model_file": "template_r9.4.1_450bps_sup.jsn", 
            "moves_out": "0", 
            "nested_output_folder": "0", 
            "noisiest_section_scaling_max_size": "8000", 
            "num_alignment_threads": "4", 
            "num_barcode_threads": "4", 
            "num_barcoding_buffers": "24", 
            "num_base_mod_threads": "2", 
            "num_callers": "4", 
            "num_extra_bases_trim": "0", 
            "num_mid_barcoding_buffers": "96", 
            "num_read_splitting_buffers": "16", 
            "num_read_splitting_threads": "4", 
            "num_reads_per_barcoding_buffer": "4", 
            "open_gap1": "40", 
            "open_gap2": "160", 
            "overlap": "100", 
            "override_scaling": "0", 
            "ping_segment_duration": "60", 
            "ping_url": "https://ping.oxfordnanoportal.com/basecall", 
            "post_out": "0", 
            "print_workflows": "0", 
            "progress_stats_frequency": "-1.000000", 
            "pt_median_offset": "2.500000", 
            "pt_minimum_read_start_index": "30", 
            "pt_required_adapter_drop": "30.000000", 
            "pt_scaling": "0", 
            "qscore_offset": "0.349800", 
            "qscore_scale": "0.972200", 
            "quiet": "0", 
            "read_batch_size": "4000", 
            "read_id_list": "", 
            "read_splitting_arrangement_files": "", 
            "read_splitting_score_matrix_filename": "", 
            "rear_window_size": "150", 
            "records_per_fastq": "0", 
            "recursive": "0", 
            "remora_models": "", 
            "require_barcodes_both_ends": "0", 
            "resume": "0", 
            "reverse_sequence": "0", 
            "sample_sheet": "", 
            "scaling_mad": "1.000000", 
            "scaling_med": "0.000000", 
            "score_matrix_filename": "5x5_mismatch_matrix.txt", 
            "start_gap1": "40", 
            "start_gap2": "40", 
            "stay_penalty": "1.000000", 
            "temp_bias": "1.000000", 
            "temp_weight": "1.000000", 
            "trace_categories_logs": "", 
            "trace_domains_config": "", 
            "trim_adapters": "0", 
            "trim_barcodes": "0", 
            "trim_min_events": "3", 
            "trim_primers": "0", 
            "trim_strategy": "dna", 
            "trim_threshold": "2.500000", 
            "u_substitution": "0", 
            "verbose_logs": "0"
        }, 
        "read_count": 5, 
        "reads_per_channel_dist": [
            {
                "channel": 40, 
                "count": 1
            }, 
            {
                "channel": 42, 
                "count": 1
            }, 
            {
                "channel": 64, 
                "count": 1
            }, 
            {
                "channel": 83, 
                "count": 1
            }, 
            {
                "channel": 86, 
                "count": 1
            }
        ], 
        "run_id": "6a94eb36563a9975b38e81c41534f6590c24867c", 
        "segment_duration": 60, 
        "segment_number": 8, 
        "segment_type": "guppy-acquisition", 
        "software": {
            "analysis": "1d_basecalling", 
            "name": "guppy-basecalling", 
            "version": "6.1.2+e0556ff"
        }, 
        "tracking_id": {
            "asic_id": "3506339664", 
            "asic_id_eeprom": "5924305", 
            "asic_temp": "35.086052", 
            "asic_version": "IA02D", 
            "auto_update": "0", 
            "auto_update_source": "https://mirror.oxfordnanoportal.com/software/MinKNOW/", 
            "bream_is_standard": "0", 
            "configuration_version": "4.2.11", 
            "device_id": "MN34875", 
            "device_type": "minion", 
            "distribution_status": "stable", 
            "distribution_version": "21.02.1", 
            "exp_script_name": "sequencing_MIN106_DNA:FLO-FLG001:SQK-LSK110", 
            "exp_script_purpose": "sequencing_run", 
            "exp_start_time": "2021-08-03T19:16:19Z", 
            "flongle_adapter_id": "FA-00993", 
            "flow_cell_id": "AHC239", 
            "flow_cell_product_code": "FLO-FLG001", 
            "guppy_version": "4.3.4+ecb2805", 
            "heatsink_temp": "35.003906", 
            "hostname": "mckennalab", 
            "installation_type": "nc", 
            "local_firmware_file": "1", 
            "msg_id": "f057581c-f7f8-4bb7-aed9-96fe79f59756", 
            "operating_system": "ubuntu 20.04", 
            "protocol_group_id": "FEX_120", 
            "protocol_run_id": "bdcd3859-b6c1-400c-b4e6-2bc287733f75", 
            "protocol_start_time": "", 
            "protocols_version": "6.1.10", 
            "run_id": "6a94eb36563a9975b38e81c41534f6590c24867c", 
            "sample_id": "96_v2", 
            "time_stamp": "2022-06-11T02:04:21Z", 
            "usb_config": "MinION_fx3_1.1.1_ONT#MinION_fpga_1.1.0#bulk#Auto", 
            "version": "4.2.5"
        }
    }, 
    {
        "aggregation": "segment", 
        "analysis_id": "e2692740-f853-4dfa-99c2-3309f56cdfdb", 
        "basecall_1d": {
            "exit_status_dist": {
                "pass": 1
            }, 
            "qscore_dist_temp": [
                {
                    "count": 1, 
                    "mean_qscore": 13.0
                }
            ], 
            "qscore_sum_temp": {
                "count": 1, 
                "mean": 13.285237312316895, 
                "sum": 13.285237312316895
            }, 
            "read_len_events_sum_temp": 26821, 
            "seq_len_bases_dist_temp": [
                {
                    "count": 1, 
                    "length": 0.0
                }
            ], 
            "seq_len_bases_sum_temp": 1, 
            "seq_len_events_dist_temp": [
                {
                    "count": 1, 
                    "length": 26000.0
                }
            ], 
            "speed_bases_per_second_dist_temp": [
                {
                    "count": 1, 
                    "speed": 1.0
                }
            ], 
            "strand_median_pa": {
                "count": 1, 
                "mean": 65.95909118652344, 
                "sum": 65.95909118652344
            }, 
            "strand_sd_pa": {
                "count": 1, 
                "mean": 8.841700553894043, 
                "sum": 8.841700553894043
            }
        }, 
        "channel_count": 1, 
        "context_tags": {
            "barcoding_enabled": "0", 
            "experiment_duration_set": "1680", 
            "experiment_type": "genomic_dna", 
            "local_basecalling": "0", 
            "package": "bream4", 
            "package_version": "6.1.10", 
            "sample_frequency": "4000", 
            "sequencing_kit": "sqk-lsk110"
        }, 
        "latest_run_time": 32395.603515625, 
        "levels_sums": {
            "count": 1, 
            "mean": 189.05520629882812, 
            "open_pore_level_sum": 189.05520629882812
        }, 
        "opts": {
            "adapter_pt_range_scale": "5.200000", 
            "additional_lamp_context_bases": "2", 
            "align_ref": "", 
            "align_type": "auto", 
            "allow_inferior_barcodes": "0", 
            "as_cpu_threads_per_scaler": "2", 
            "as_gpu_runners_per_device": "2", 
            "as_model_file": "", 
            "as_num_scalers": "4", 
            "as_reads_per_runner": "32", 
            "bam_methylation_threshold": "5.000000", 
            "bam_out": "0", 
            "barcode_kits": "", 
            "barcode_nested_output_folder": "0", 
            "beam_cut": "100.000000", 
            "beam_width": "32", 
            "bed_file": "", 
            "builtin_scripts": "1", 
            "calib_detect": "0", 
            "calib_max_sequence_length": "3800", 
            "calib_min_coverage": "0.600000", 
            "calib_min_sequence_length": "3000", 
            "calib_reference": "lambda_3.6kb.fasta", 
            "chunk_size": "2000", 
            "chunks_per_caller": "10000", 
            "chunks_per_runner": "208", 
            "client_id": "-1", 
            "compress_fastq": "1", 
            "cpu_threads_per_caller": "4", 
            "detect_adapter": "0", 
            "detect_barcodes": "0", 
            "detect_mid_strand_adapter": "0", 
            "detect_mid_strand_barcodes": "0", 
            "detect_primer": "0", 
            "device": "cuda:0", 
            "disable_pings": "0", 
            "disable_qscore_filtering": "0", 
            "dmean_threshold": "1.000000", 
            "dmean_win_size": "2", 
            "do_read_splitting": "0", 
            "duplex_window_size_max": "1000", 
            "duplex_window_size_min": "200", 
            "end_gap1": "40", 
            "end_gap2": "40", 
            "extend_gap1": "40", 
            "extend_gap2": "160", 
            "fast5_out": "0", 
            "flowcell": "", 
            "front_window_size": "150", 
            "gpu_runners_per_device": "12", 
            "high_priority_threshold": "10", 
            "index": "0", 
            "input_file_list": "", 
            "int8_mode": "0", 
            "jump_threshold": "1.000000", 
            "kernel_path": "", 
            "kit": "", 
            "lamp_kit": "", 
            "log_speed_frequency": "0", 
            "max_queued_reads": "2000", 
            "max_read_split_depth": "2", 
            "max_search_len": "1000", 
            "medium_priority_threshold": "4", 
            "min_length_lamp_context": "30", 
            "min_length_lamp_target": "70", 
            "min_qscore": "10.000000", 
            "min_score_adapter": "60.000000", 
            "min_score_adapter_mid": "50.000000", 
            "min_score_barcode_front": "60.000000", 
            "min_score_barcode_mask": "40.000000", 
            "min_score_barcode_mid": "50.000000", 
            "min_score_barcode_rear": "60.000000", 
            "min_score_lamp": "80.000000", 
            "min_score_lamp_mask": "50.000000", 
            "min_score_lamp_target": "50.000000", 
            "min_score_primer": "60.000000", 
            "min_score_read_splitting": "70.000000", 
            "model_file": "template_r9.4.1_450bps_sup.jsn", 
            "moves_out": "0", 
            "nested_output_folder": "0", 
            "noisiest_section_scaling_max_size": "8000", 
            "num_alignment_threads": "4", 
            "num_barcode_threads": "4", 
            "num_barcoding_buffers": "24", 
            "num_base_mod_threads": "2", 
            "num_callers": "4", 
            "num_extra_bases_trim": "0", 
            "num_mid_barcoding_buffers": "96", 
            "num_read_splitting_buffers": "16", 
            "num_read_splitting_threads": "4", 
            "num_reads_per_barcoding_buffer": "4", 
            "open_gap1": "40", 
            "open_gap2": "160", 
            "overlap": "100", 
            "override_scaling": "0", 
            "ping_segment_duration": "60", 
            "ping_url": "https://ping.oxfordnanoportal.com/basecall", 
            "post_out": "0", 
            "print_workflows": "0", 
            "progress_stats_frequency": "-1.000000", 
            "pt_median_offset": "2.500000", 
            "pt_minimum_read_start_index": "30", 
            "pt_required_adapter_drop": "30.000000", 
            "pt_scaling": "0", 
            "qscore_offset": "0.349800", 
            "qscore_scale": "0.972200", 
            "quiet": "0", 
            "read_batch_size": "4000", 
            "read_id_list": "", 
            "read_splitting_arrangement_files": "", 
            "read_splitting_score_matrix_filename": "", 
            "rear_window_size": "150", 
            "records_per_fastq": "0", 
            "recursive": "0", 
            "remora_models": "", 
            "require_barcodes_both_ends": "0", 
            "resume": "0", 
            "reverse_sequence": "0", 
            "sample_sheet": "", 
            "scaling_mad": "1.000000", 
            "scaling_med": "0.000000", 
            "score_matrix_filename": "5x5_mismatch_matrix.txt", 
            "start_gap1": "40", 
            "start_gap2": "40", 
            "stay_penalty": "1.000000", 
            "temp_bias": "1.000000", 
            "temp_weight": "1.000000", 
            "trace_categories_logs": "", 
            "trace_domains_config": "", 
            "trim_adapters": "0", 
            "trim_barcodes": "0", 
            "trim_min_events": "3", 
            "trim_primers": "0", 
            "trim_strategy": "dna", 
            "trim_threshold": "2.500000", 
            "u_substitution": "0", 
            "verbose_logs": "0"
        }, 
        "read_count": 1, 
        "reads_per_channel_dist": [
            {
                "channel": 89, 
                "count": 1
            }
        ], 
        "run_id": "6a94eb36563a9975b38e81c41534f6590c24867c", 
        "segment_duration": 60, 
        "segment_number": 9, 
        "segment_type": "guppy-acquisition", 
        "software": {
            "analysis": "1d_basecalling", 
            "name": "guppy-basecalling", 
            "version": "6.1.2+e0556ff"
        }, 
        "tracking_id": {
            "asic_id": "3506339664", 
            "asic_id_eeprom": "5924305", 
            "asic_temp": "35.086052", 
            "asic_version": "IA02D", 
            "auto_update": "0", 
            "auto_update_source": "https://mirror.oxfordnanoportal.com/software/MinKNOW/", 
            "bream_is_standard": "0", 
            "configuration_version": "4.2.11", 
            "device_id": "MN34875", 
            "device_type": "minion", 
            "distribution_status": "stable", 
            "distribution_version": "21.02.1", 
            "exp_script_name": "sequencing_MIN106_DNA:FLO-FLG001:SQK-LSK110", 
            "exp_script_purpose": "sequencing_run", 
            "exp_start_time": "2021-08-03T19:16:19Z", 
            "flongle_adapter_id": "FA-00993", 
            "flow_cell_id": "AHC239", 
            "flow_cell_product_code": "FLO-FLG001", 
            "guppy_version": "4.3.4+ecb2805", 
            "heatsink_temp": "35.003906", 
            "hostname": "mckennalab", 
            "installation_type": "nc", 
            "local_firmware_file": "1", 
            "msg_id": "9724cc36-b2ca-4250-bd76-d294daadff20", 
            "operating_system": "ubuntu 20.04", 
            "protocol_group_id": "FEX_120", 
            "protocol_run_id": "bdcd3859-b6c1-400c-b4e6-2bc287733f75", 
            "protocol_start_time": "", 
            "protocols_version": "6.1.10", 
            "run_id": "6a94eb36563a9975b38e81c41534f6590c24867c", 
            "sample_id": "96_v2", 
            "time_stamp": "2022-06-11T02:04:21Z", 
            "usb_config": "MinION_fx3_1.1.1_ONT#MinION_fpga_1.1.0#bulk#Auto", 
            "version": "4.2.5"
        }
    }, 
    {
        "aggregation": "segment", 
        "analysis_id": "e2692740-f853-4dfa-99c2-3309f56cdfdb", 
        "basecall_1d": {
            "exit_status_dist": {
                "pass": 1
            }, 
            "qscore_dist_temp": [
                {
                    "count": 1, 
                    "mean_qscore": 13.5
                }
            ], 
            "qscore_sum_temp": {
                "count": 1, 
                "mean": 13.805153846740723, 
                "sum": 13.805153846740723
            }, 
            "read_len_events_sum_temp": 23415, 
            "seq_len_bases_dist_temp": [
                {
                    "count": 1, 
                    "length": 0.0
                }
            ], 
            "seq_len_bases_sum_temp": 1, 
            "seq_len_events_dist_temp": [
                {
                    "count": 1, 
                    "length": 23000.0
                }
            ], 
            "speed_bases_per_second_dist_temp": [
                {
                    "count": 1, 
                    "speed": 1.0
                }
            ], 
            "strand_median_pa": {
                "count": 1, 
                "mean": 67.72743225097656, 
                "sum": 67.72743225097656
            }, 
            "strand_sd_pa": {
                "count": 1, 
                "mean": 9.018534660339355, 
                "sum": 9.018534660339355
            }
        }, 
        "channel_count": 1, 
        "context_tags": {
            "barcoding_enabled": "0", 
            "experiment_duration_set": "1680", 
            "experiment_type": "genomic_dna", 
            "local_basecalling": "0", 
            "package": "bream4", 
            "package_version": "6.1.10", 
            "sample_frequency": "4000", 
            "sequencing_kit": "sqk-lsk110"
        }, 
        "latest_run_time": 34076.3203125, 
        "levels_sums": {
            "count": 1, 
            "mean": 187.48468017578125, 
            "open_pore_level_sum": 187.48468017578125
        }, 
        "opts": {
            "adapter_pt_range_scale": "5.200000", 
            "additional_lamp_context_bases": "2", 
            "align_ref": "", 
            "align_type": "auto", 
            "allow_inferior_barcodes": "0", 
            "as_cpu_threads_per_scaler": "2", 
            "as_gpu_runners_per_device": "2", 
            "as_model_file": "", 
            "as_num_scalers": "4", 
            "as_reads_per_runner": "32", 
            "bam_methylation_threshold": "5.000000", 
            "bam_out": "0", 
            "barcode_kits": "", 
            "barcode_nested_output_folder": "0", 
            "beam_cut": "100.000000", 
            "beam_width": "32", 
            "bed_file": "", 
            "builtin_scripts": "1", 
            "calib_detect": "0", 
            "calib_max_sequence_length": "3800", 
            "calib_min_coverage": "0.600000", 
            "calib_min_sequence_length": "3000", 
            "calib_reference": "lambda_3.6kb.fasta", 
            "chunk_size": "2000", 
            "chunks_per_caller": "10000", 
            "chunks_per_runner": "208", 
            "client_id": "-1", 
            "compress_fastq": "1", 
            "cpu_threads_per_caller": "4", 
            "detect_adapter": "0", 
            "detect_barcodes": "0", 
            "detect_mid_strand_adapter": "0", 
            "detect_mid_strand_barcodes": "0", 
            "detect_primer": "0", 
            "device": "cuda:0", 
            "disable_pings": "0", 
            "disable_qscore_filtering": "0", 
            "dmean_threshold": "1.000000", 
            "dmean_win_size": "2", 
            "do_read_splitting": "0", 
            "duplex_window_size_max": "1000", 
            "duplex_window_size_min": "200", 
            "end_gap1": "40", 
            "end_gap2": "40", 
            "extend_gap1": "40", 
            "extend_gap2": "160", 
            "fast5_out": "0", 
            "flowcell": "", 
            "front_window_size": "150", 
            "gpu_runners_per_device": "12", 
            "high_priority_threshold": "10", 
            "index": "0", 
            "input_file_list": "", 
            "int8_mode": "0", 
            "jump_threshold": "1.000000", 
            "kernel_path": "", 
            "kit": "", 
            "lamp_kit": "", 
            "log_speed_frequency": "0", 
            "max_queued_reads": "2000", 
            "max_read_split_depth": "2", 
            "max_search_len": "1000", 
            "medium_priority_threshold": "4", 
            "min_length_lamp_context": "30", 
            "min_length_lamp_target": "70", 
            "min_qscore": "10.000000", 
            "min_score_adapter": "60.000000", 
            "min_score_adapter_mid": "50.000000", 
            "min_score_barcode_front": "60.000000", 
            "min_score_barcode_mask": "40.000000", 
            "min_score_barcode_mid": "50.000000", 
            "min_score_barcode_rear": "60.000000", 
            "min_score_lamp": "80.000000", 
            "min_score_lamp_mask": "50.000000", 
            "min_score_lamp_target": "50.000000", 
            "min_score_primer": "60.000000", 
            "min_score_read_splitting": "70.000000", 
            "model_file": "template_r9.4.1_450bps_sup.jsn", 
            "moves_out": "0", 
            "nested_output_folder": "0", 
            "noisiest_section_scaling_max_size": "8000", 
            "num_alignment_threads": "4", 
            "num_barcode_threads": "4", 
            "num_barcoding_buffers": "24", 
            "num_base_mod_threads": "2", 
            "num_callers": "4", 
            "num_extra_bases_trim": "0", 
            "num_mid_barcoding_buffers": "96", 
            "num_read_splitting_buffers": "16", 
            "num_read_splitting_threads": "4", 
            "num_reads_per_barcoding_buffer": "4", 
            "open_gap1": "40", 
            "open_gap2": "160", 
            "overlap": "100", 
            "override_scaling": "0", 
            "ping_segment_duration": "60", 
            "ping_url": "https://ping.oxfordnanoportal.com/basecall", 
            "post_out": "0", 
            "print_workflows": "0", 
            "progress_stats_frequency": "-1.000000", 
            "pt_median_offset": "2.500000", 
            "pt_minimum_read_start_index": "30", 
            "pt_required_adapter_drop": "30.000000", 
            "pt_scaling": "0", 
            "qscore_offset": "0.349800", 
            "qscore_scale": "0.972200", 
            "quiet": "0", 
            "read_batch_size": "4000", 
            "read_id_list": "", 
            "read_splitting_arrangement_files": "", 
            "read_splitting_score_matrix_filename": "", 
            "rear_window_size": "150", 
            "records_per_fastq": "0", 
            "recursive": "0", 
            "remora_models": "", 
            "require_barcodes_both_ends": "0", 
            "resume": "0", 
            "reverse_sequence": "0", 
            "sample_sheet": "", 
            "scaling_mad": "1.000000", 
            "scaling_med": "0.000000", 
            "score_matrix_filename": "5x5_mismatch_matrix.txt", 
            "start_gap1": "40", 
            "start_gap2": "40", 
            "stay_penalty": "1.000000", 
            "temp_bias": "1.000000", 
            "temp_weight": "1.000000", 
            "trace_categories_logs": "", 
            "trace_domains_config": "", 
            "trim_adapters": "0", 
            "trim_barcodes": "0", 
            "trim_min_events": "3", 
            "trim_primers": "0", 
            "trim_strategy": "dna", 
            "trim_threshold": "2.500000", 
            "u_substitution": "0", 
            "verbose_logs": "0"
        }, 
        "read_count": 1, 
        "reads_per_channel_dist": [
            {
                "channel": 49, 
                "count": 1
            }
        ], 
        "run_id": "6a94eb36563a9975b38e81c41534f6590c24867c", 
        "segment_duration": 60, 
        "segment_number": 10, 
        "segment_type": "guppy-acquisition", 
        "software": {
            "analysis": "1d_basecalling", 
            "name": "guppy-basecalling", 
            "version": "6.1.2+e0556ff"
        }, 
        "tracking_id": {
            "asic_id": "3506339664", 
            "asic_id_eeprom": "5924305", 
            "asic_temp": "35.086052", 
            "asic_version": "IA02D", 
            "auto_update": "0", 
            "auto_update_source": "https://mirror.oxfordnanoportal.com/software/MinKNOW/", 
            "bream_is_standard": "0", 
            "configuration_version": "4.2.11", 
            "device_id": "MN34875", 
            "device_type": "minion", 
            "distribution_status": "stable", 
            "distribution_version": "21.02.1", 
            "exp_script_name": "sequencing_MIN106_DNA:FLO-FLG001:SQK-LSK110", 
            "exp_script_purpose": "sequencing_run", 
            "exp_start_time": "2021-08-03T19:16:19Z", 
            "flongle_adapter_id": "FA-00993", 
            "flow_cell_id": "AHC239", 
            "flow_cell_product_code": "FLO-FLG001", 
            "guppy_version": "4.3.4+ecb2805", 
            "heatsink_temp": "35.003906", 
            "hostname": "mckennalab", 
            "installation_type": "nc", 
            "local_firmware_file": "1", 
            "msg_id": "f65b71c1-3431-49dd-9a13-1ad2294f0be9", 
            "operating_system": "ubuntu 20.04", 
            "protocol_group_id": "FEX_120", 
            "protocol_run_id": "bdcd3859-b6c1-400c-b4e6-2bc287733f75", 
            "protocol_start_time": "", 
            "protocols_version": "6.1.10", 
            "run_id": "6a94eb36563a9975b38e81c41534f6590c24867c", 
            "sample_id": "96_v2", 
            "time_stamp": "2022-06-11T02:04:21Z", 
            "usb_config": "MinION_fx3_1.1.1_ONT#MinION_fpga_1.1.0#bulk#Auto", 
            "version": "4.2.5"
        }
    }, 
    {
        "aggregation": "segment", 
        "analysis_id": "e2692740-f853-4dfa-99c2-3309f56cdfdb", 
        "basecall_1d": {
            "exit_status_dist": {
                "pass": 1
            }, 
            "qscore_dist_temp": [
                {
                    "count": 1, 
                    "mean_qscore": 10.5
                }
            ], 
            "qscore_sum_temp": {
                "count": 1, 
                "mean": 10.806424140930176, 
                "sum": 10.806424140930176
            }, 
            "read_len_events_sum_temp": 20434, 
            "seq_len_bases_dist_temp": [
                {
                    "count": 1, 
                    "length": 0.0
                }
            ], 
            "seq_len_bases_sum_temp": 1, 
            "seq_len_events_dist_temp": [
                {
                    "count": 1, 
                    "length": 20000.0
                }
            ], 
            "speed_bases_per_second_dist_temp": [
                {
                    "count": 1, 
                    "speed": 1.0
                }
            ], 
            "strand_median_pa": {
                "count": 1, 
                "mean": 70.2031021118164, 
                "sum": 70.2031021118164
            }, 
            "strand_sd_pa": {
                "count": 1, 
                "mean": 9.018534660339355, 
                "sum": 9.018534660339355
            }
        }, 
        "channel_count": 1, 
        "context_tags": {
            "barcoding_enabled": "0", 
            "experiment_duration_set": "1680", 
            "experiment_type": "genomic_dna", 
            "local_basecalling": "0", 
            "package": "bream4", 
            "package_version": "6.1.10", 
            "sample_frequency": "4000", 
            "sequencing_kit": "sqk-lsk110"
        }, 
        "latest_run_time": 36270.02734375, 
        "levels_sums": {
            "count": 1, 
            "mean": 195.75526428222656, 
            "open_pore_level_sum": 195.75526428222656
        }, 
        "opts": {
            "adapter_pt_range_scale": "5.200000", 
            "additional_lamp_context_bases": "2", 
            "align_ref": "", 
            "align_type": "auto", 
            "allow_inferior_barcodes": "0", 
            "as_cpu_threads_per_scaler": "2", 
            "as_gpu_runners_per_device": "2", 
            "as_model_file": "", 
            "as_num_scalers": "4", 
            "as_reads_per_runner": "32", 
            "bam_methylation_threshold": "5.000000", 
            "bam_out": "0", 
            "barcode_kits": "", 
            "barcode_nested_output_folder": "0", 
            "beam_cut": "100.000000", 
            "beam_width": "32", 
            "bed_file": "", 
            "builtin_scripts": "1", 
            "calib_detect": "0", 
            "calib_max_sequence_length": "3800", 
            "calib_min_coverage": "0.600000", 
            "calib_min_sequence_length": "3000", 
            "calib_reference": "lambda_3.6kb.fasta", 
            "chunk_size": "2000", 
            "chunks_per_caller": "10000", 
            "chunks_per_runner": "208", 
            "client_id": "-1", 
            "compress_fastq": "1", 
            "cpu_threads_per_caller": "4", 
            "detect_adapter": "0", 
            "detect_barcodes": "0", 
            "detect_mid_strand_adapter": "0", 
            "detect_mid_strand_barcodes": "0", 
            "detect_primer": "0", 
            "device": "cuda:0", 
            "disable_pings": "0", 
            "disable_qscore_filtering": "0", 
            "dmean_threshold": "1.000000", 
            "dmean_win_size": "2", 
            "do_read_splitting": "0", 
            "duplex_window_size_max": "1000", 
            "duplex_window_size_min": "200", 
            "end_gap1": "40", 
            "end_gap2": "40", 
            "extend_gap1": "40", 
            "extend_gap2": "160", 
            "fast5_out": "0", 
            "flowcell": "", 
            "front_window_size": "150", 
            "gpu_runners_per_device": "12", 
            "high_priority_threshold": "10", 
            "index": "0", 
            "input_file_list": "", 
            "int8_mode": "0", 
            "jump_threshold": "1.000000", 
            "kernel_path": "", 
            "kit": "", 
            "lamp_kit": "", 
            "log_speed_frequency": "0", 
            "max_queued_reads": "2000", 
            "max_read_split_depth": "2", 
            "max_search_len": "1000", 
            "medium_priority_threshold": "4", 
            "min_length_lamp_context": "30", 
            "min_length_lamp_target": "70", 
            "min_qscore": "10.000000", 
            "min_score_adapter": "60.000000", 
            "min_score_adapter_mid": "50.000000", 
            "min_score_barcode_front": "60.000000", 
            "min_score_barcode_mask": "40.000000", 
            "min_score_barcode_mid": "50.000000", 
            "min_score_barcode_rear": "60.000000", 
            "min_score_lamp": "80.000000", 
            "min_score_lamp_mask": "50.000000", 
            "min_score_lamp_target": "50.000000", 
            "min_score_primer": "60.000000", 
            "min_score_read_splitting": "70.000000", 
            "model_file": "template_r9.4.1_450bps_sup.jsn", 
            "moves_out": "0", 
            "nested_output_folder": "0", 
            "noisiest_section_scaling_max_size": "8000", 
            "num_alignment_threads": "4", 
            "num_barcode_threads": "4", 
            "num_barcoding_buffers": "24", 
            "num_base_mod_threads": "2", 
            "num_callers": "4", 
            "num_extra_bases_trim": "0", 
            "num_mid_barcoding_buffers": "96", 
            "num_read_splitting_buffers": "16", 
            "num_read_splitting_threads": "4", 
            "num_reads_per_barcoding_buffer": "4", 
            "open_gap1": "40", 
            "open_gap2": "160", 
            "overlap": "100", 
            "override_scaling": "0", 
            "ping_segment_duration": "60", 
            "ping_url": "https://ping.oxfordnanoportal.com/basecall", 
            "post_out": "0", 
            "print_workflows": "0", 
            "progress_stats_frequency": "-1.000000", 
            "pt_median_offset": "2.500000", 
            "pt_minimum_read_start_index": "30", 
            "pt_required_adapter_drop": "30.000000", 
            "pt_scaling": "0", 
            "qscore_offset": "0.349800", 
            "qscore_scale": "0.972200", 
            "quiet": "0", 
            "read_batch_size": "4000", 
            "read_id_list": "", 
            "read_splitting_arrangement_files": "", 
            "read_splitting_score_matrix_filename": "", 
            "rear_window_size": "150", 
            "records_per_fastq": "0", 
            "recursive": "0", 
            "remora_models": "", 
            "require_barcodes_both_ends": "0", 
            "resume": "0", 
            "reverse_sequence": "0", 
            "sample_sheet": "", 
            "scaling_mad": "1.000000", 
            "scaling_med": "0.000000", 
            "score_matrix_filename": "5x5_mismatch_matrix.txt", 
            "start_gap1": "40", 
            "start_gap2": "40", 
            "stay_penalty": "1.000000", 
            "temp_bias": "1.000000", 
            "temp_weight": "1.000000", 
            "trace_categories_logs": "", 
            "trace_domains_config": "", 
            "trim_adapters": "0", 
            "trim_barcodes": "0", 
            "trim_min_events": "3", 
            "trim_primers": "0", 
            "trim_strategy": "dna", 
            "trim_threshold": "2.500000", 
            "u_substitution": "0", 
            "verbose_logs": "0"
        }, 
        "read_count": 1, 
        "reads_per_channel_dist": [
            {
                "channel": 62, 
                "count": 1
            }
        ], 
        "run_id": "6a94eb36563a9975b38e81c41534f6590c24867c", 
        "segment_duration": 60, 
        "segment_number": 11, 
        "segment_type": "guppy-acquisition", 
        "software": {
            "analysis": "1d_basecalling", 
            "name": "guppy-basecalling", 
            "version": "6.1.2+e0556ff"
        }, 
        "tracking_id": {
            "asic_id": "3506339664", 
            "asic_id_eeprom": "5924305", 
            "asic_temp": "35.086052", 
            "asic_version": "IA02D", 
            "auto_update": "0", 
            "auto_update_source": "https://mirror.oxfordnanoportal.com/software/MinKNOW/", 
            "bream_is_standard": "0", 
            "configuration_version": "4.2.11", 
            "device_id": "MN34875", 
            "device_type": "minion", 
            "distribution_status": "stable", 
            "distribution_version": "21.02.1", 
            "exp_script_name": "sequencing_MIN106_DNA:FLO-FLG001:SQK-LSK110", 
            "exp_script_purpose": "sequencing_run", 
            "exp_start_time": "2021-08-03T19:16:19Z", 
            "flongle_adapter_id": "FA-00993", 
            "flow_cell_id": "AHC239", 
            "flow_cell_product_code": "FLO-FLG001", 
            "guppy_version": "4.3.4+ecb2805", 
            "heatsink_temp": "35.003906", 
            "hostname": "mckennalab", 
            "installation_type": "nc", 
            "local_firmware_file": "1", 
            "msg_id": "03d3f493-3db6-4949-a591-9b7113ebced4", 
            "operating_system": "ubuntu 20.04", 
            "protocol_group_id": "FEX_120", 
            "protocol_run_id": "bdcd3859-b6c1-400c-b4e6-2bc287733f75", 
            "protocol_start_time": "", 
            "protocols_version": "6.1.10", 
            "run_id": "6a94eb36563a9975b38e81c41534f6590c24867c", 
            "sample_id": "96_v2", 
            "time_stamp": "2022-06-11T02:04:21Z", 
            "usb_config": "MinION_fx3_1.1.1_ONT#MinION_fpga_1.1.0#bulk#Auto", 
            "version": "4.2.5"
        }
    }, 
    {
        "aggregation": "segment", 
        "analysis_id": "e2692740-f853-4dfa-99c2-3309f56cdfdb", 
        "basecall_1d": {
            "exit_status_dist": {
                "pass": 4
            }, 
            "qscore_dist_temp": [
                {
                    "count": 1, 
                    "mean_qscore": 11.0
                }, 
                {
                    "count": 1, 
                    "mean_qscore": 12.5
                }, 
                {
                    "count": 1, 
                    "mean_qscore": 13.0
                }, 
                {
                    "count": 1, 
                    "mean_qscore": 14.5
                }
            ], 
            "qscore_sum_temp": {
                "count": 4, 
                "mean": 12.896550178527832, 
                "sum": 51.58620071411133
            }, 
            "read_len_events_sum_temp": 119180, 
            "seq_len_bases_dist_temp": [
                {
                    "count": 4, 
                    "length": 0.0
                }
            ], 
            "seq_len_bases_sum_temp": 4, 
            "seq_len_events_dist_temp": [
                {
                    "count": 2, 
                    "length": 26000.0
                }, 
                {
                    "count": 1, 
                    "length": 27000.0
                }, 
                {
                    "count": 1, 
                    "length": 38000.0
                }
            ], 
            "speed_bases_per_second_dist_temp": [
                {
                    "count": 4, 
                    "speed": 1.0
                }
            ], 
            "strand_median_pa": {
                "count": 4, 
                "mean": 70.822021484375, 
                "sum": 283.2880859375
            }, 
            "strand_sd_pa": {
                "count": 4, 
                "mean": 9.327994346618652, 
                "sum": 37.31197738647461
            }
        }, 
        "channel_count": 4, 
        "context_tags": {
            "barcoding_enabled": "0", 
            "experiment_duration_set": "1680", 
            "experiment_type": "genomic_dna", 
            "local_basecalling": "0", 
            "package": "bream4", 
            "package_version": "6.1.10", 
            "sample_frequency": "4000", 
            "sequencing_kit": "sqk-lsk110"
        }, 
        "latest_run_time": 41668.8671875, 
        "levels_sums": {
            "count": 4, 
            "mean": 198.3684539794922, 
            "open_pore_level_sum": 793.47381591796875
        }, 
        "opts": {
            "adapter_pt_range_scale": "5.200000", 
            "additional_lamp_context_bases": "2", 
            "align_ref": "", 
            "align_type": "auto", 
            "allow_inferior_barcodes": "0", 
            "as_cpu_threads_per_scaler": "2", 
            "as_gpu_runners_per_device": "2", 
            "as_model_file": "", 
            "as_num_scalers": "4", 
            "as_reads_per_runner": "32", 
            "bam_methylation_threshold": "5.000000", 
            "bam_out": "0", 
            "barcode_kits": "", 
            "barcode_nested_output_folder": "0", 
            "beam_cut": "100.000000", 
            "beam_width": "32", 
            "bed_file": "", 
            "builtin_scripts": "1", 
            "calib_detect": "0", 
            "calib_max_sequence_length": "3800", 
            "calib_min_coverage": "0.600000", 
            "calib_min_sequence_length": "3000", 
            "calib_reference": "lambda_3.6kb.fasta", 
            "chunk_size": "2000", 
            "chunks_per_caller": "10000", 
            "chunks_per_runner": "208", 
            "client_id": "-1", 
            "compress_fastq": "1", 
            "cpu_threads_per_caller": "4", 
            "detect_adapter": "0", 
            "detect_barcodes": "0", 
            "detect_mid_strand_adapter": "0", 
            "detect_mid_strand_barcodes": "0", 
            "detect_primer": "0", 
            "device": "cuda:0", 
            "disable_pings": "0", 
            "disable_qscore_filtering": "0", 
            "dmean_threshold": "1.000000", 
            "dmean_win_size": "2", 
            "do_read_splitting": "0", 
            "duplex_window_size_max": "1000", 
            "duplex_window_size_min": "200", 
            "end_gap1": "40", 
            "end_gap2": "40", 
            "extend_gap1": "40", 
            "extend_gap2": "160", 
            "fast5_out": "0", 
            "flowcell": "", 
            "front_window_size": "150", 
            "gpu_runners_per_device": "12", 
            "high_priority_threshold": "10", 
            "index": "0", 
            "input_file_list": "", 
            "int8_mode": "0", 
            "jump_threshold": "1.000000", 
            "kernel_path": "", 
            "kit": "", 
            "lamp_kit": "", 
            "log_speed_frequency": "0", 
            "max_queued_reads": "2000", 
            "max_read_split_depth": "2", 
            "max_search_len": "1000", 
            "medium_priority_threshold": "4", 
            "min_length_lamp_context": "30", 
            "min_length_lamp_target": "70", 
            "min_qscore": "10.000000", 
            "min_score_adapter": "60.000000", 
            "min_score_adapter_mid": "50.000000", 
            "min_score_barcode_front": "60.000000", 
            "min_score_barcode_mask": "40.000000", 
            "min_score_barcode_mid": "50.000000", 
            "min_score_barcode_rear": "60.000000", 
            "min_score_lamp": "80.000000", 
            "min_score_lamp_mask": "50.000000", 
            "min_score_lamp_target": "50.000000", 
            "min_score_primer": "60.000000", 
            "min_score_read_splitting": "70.000000", 
            "model_file": "template_r9.4.1_450bps_sup.jsn", 
            "moves_out": "0", 
            "nested_output_folder": "0", 
            "noisiest_section_scaling_max_size": "8000", 
            "num_alignment_threads": "4", 
            "num_barcode_threads": "4", 
            "num_barcoding_buffers": "24", 
            "num_base_mod_threads": "2", 
            "num_callers": "4", 
            "num_extra_bases_trim": "0", 
            "num_mid_barcoding_buffers": "96", 
            "num_read_splitting_buffers": "16", 
            "num_read_splitting_threads": "4", 
            "num_reads_per_barcoding_buffer": "4", 
            "open_gap1": "40", 
            "open_gap2": "160", 
            "overlap": "100", 
            "override_scaling": "0", 
            "ping_segment_duration": "60", 
            "ping_url": "https://ping.oxfordnanoportal.com/basecall", 
            "post_out": "0", 
            "print_workflows": "0", 
            "progress_stats_frequency": "-1.000000", 
            "pt_median_offset": "2.500000", 
            "pt_minimum_read_start_index": "30", 
            "pt_required_adapter_drop": "30.000000", 
            "pt_scaling": "0", 
            "qscore_offset": "0.349800", 
            "qscore_scale": "0.972200", 
            "quiet": "0", 
            "read_batch_size": "4000", 
            "read_id_list": "", 
            "read_splitting_arrangement_files": "", 
            "read_splitting_score_matrix_filename": "", 
            "rear_window_size": "150", 
            "records_per_fastq": "0", 
            "recursive": "0", 
            "remora_models": "", 
            "require_barcodes_both_ends": "0", 
            "resume": "0", 
            "reverse_sequence": "0", 
            "sample_sheet": "", 
            "scaling_mad": "1.000000", 
            "scaling_med": "0.000000", 
            "score_matrix_filename": "5x5_mismatch_matrix.txt", 
            "start_gap1": "40", 
            "start_gap2": "40", 
            "stay_penalty": "1.000000", 
            "temp_bias": "1.000000", 
            "temp_weight": "1.000000", 
            "trace_categories_logs": "", 
            "trace_domains_config": "", 
            "trim_adapters": "0", 
            "trim_barcodes": "0", 
            "trim_min_events": "3", 
            "trim_primers": "0", 
            "trim_strategy": "dna", 
            "trim_threshold": "2.500000", 
            "u_substitution": "0", 
            "verbose_logs": "0"
        }, 
        "read_count": 4, 
        "reads_per_channel_dist": [
            {
                "channel": 16, 
                "count": 1
            }, 
            {
                "channel": 54, 
                "count": 1
            }, 
            {
                "channel": 64, 
                "count": 1
            }, 
            {
                "channel": 111, 
                "count": 1
            }
        ], 
        "run_id": "6a94eb36563a9975b38e81c41534f6590c24867c", 
        "segment_duration": 60, 
        "segment_number": 12, 
        "segment_type": "guppy-acquisition", 
        "software": {
            "analysis": "1d_basecalling", 
            "name": "guppy-basecalling", 
            "version": "6.1.2+e0556ff"
        }, 
        "tracking_id": {
            "asic_id": "3506339664", 
            "asic_id_eeprom": "5924305", 
            "asic_temp": "35.086052", 
            "asic_version": "IA02D", 
            "auto_update": "0", 
            "auto_update_source": "https://mirror.oxfordnanoportal.com/software/MinKNOW/", 
            "bream_is_standard": "0", 
            "configuration_version": "4.2.11", 
            "device_id": "MN34875", 
            "device_type": "minion", 
            "distribution_status": "stable", 
            "distribution_version": "21.02.1", 
            "exp_script_name": "sequencing_MIN106_DNA:FLO-FLG001:SQK-LSK110", 
            "exp_script_purpose": "sequencing_run", 
            "exp_start_time": "2021-08-03T19:16:19Z", 
            "flongle_adapter_id": "FA-00993", 
            "flow_cell_id": "AHC239", 
            "flow_cell_product_code": "FLO-FLG001", 
            "guppy_version": "4.3.4+ecb2805", 
            "heatsink_temp": "35.003906", 
            "hostname": "mckennalab", 
            "installation_type": "nc", 
            "local_firmware_file": "1", 
            "msg_id": "7a04cb4b-901b-4d99-b53f-c038f69d78b6", 
            "operating_system": "ubuntu 20.04", 
            "protocol_group_id": "FEX_120", 
            "protocol_run_id": "bdcd3859-b6c1-400c-b4e6-2bc287733f75", 
            "protocol_start_time": "", 
            "protocols_version": "6.1.10", 
            "run_id": "6a94eb36563a9975b38e81c41534f6590c24867c", 
            "sample_id": "96_v2", 
            "time_stamp": "2022-06-11T02:04:21Z", 
            "usb_config": "MinION_fx3_1.1.1_ONT#MinION_fpga_1.1.0#bulk#Auto", 
            "version": "4.2.5"
        }
    }, 
    {
        "aggregation": "segment", 
        "analysis_id": "e2692740-f853-4dfa-99c2-3309f56cdfdb", 
        "basecall_1d": {
            "exit_status_dist": {
                "pass": 2
            }, 
            "qscore_dist_temp": [
                {
                    "count": 1, 
                    "mean_qscore": 12.0
                }, 
                {
                    "count": 1, 
                    "mean_qscore": 13.5
                }
            ], 
            "qscore_sum_temp": {
                "count": 2, 
                "mean": 13.057008743286133, 
                "sum": 26.114017486572266
            }, 
            "read_len_events_sum_temp": 64136, 
            "seq_len_bases_dist_temp": [
                {
                    "count": 2, 
                    "length": 0.0
                }
            ], 
            "seq_len_bases_sum_temp": 2, 
            "seq_len_events_dist_temp": [
                {
                    "count": 1, 
                    "length": 23000.0
                }, 
                {
                    "count": 1, 
                    "length": 40000.0
                }
            ], 
            "speed_bases_per_second_dist_temp": [
                {
                    "count": 2, 
                    "speed": 1.0
                }
            ], 
            "strand_median_pa": {
                "count": 2, 
                "mean": 69.93785095214844, 
                "sum": 139.87570190429688
            }, 
            "strand_sd_pa": {
                "count": 2, 
                "mean": 9.37220287322998, 
                "sum": 18.74440574645996
            }
        }, 
        "channel_count": 2, 
        "context_tags": {
            "barcoding_enabled": "0", 
            "experiment_duration_set": "1680", 
            "experiment_type": "genomic_dna", 
            "local_basecalling": "0", 
            "package": "bream4", 
            "package_version": "6.1.10", 
            "sample_frequency": "4000", 
            "sequencing_kit": "sqk-lsk110"
        }, 
        "latest_run_time": 46200.4921875, 
        "levels_sums": {
            "count": 2, 
            "mean": 193.0042266845703, 
            "open_pore_level_sum": 386.0084533691406
        }, 
        "opts": {
            "adapter_pt_range_scale": "5.200000", 
            "additional_lamp_context_bases": "2", 
            "align_ref": "", 
            "align_type": "auto", 
            "allow_inferior_barcodes": "0", 
            "as_cpu_threads_per_scaler": "2", 
            "as_gpu_runners_per_device": "2", 
            "as_model_file": "", 
            "as_num_scalers": "4", 
            "as_reads_per_runner": "32", 
            "bam_methylation_threshold": "5.000000", 
            "bam_out": "0", 
            "barcode_kits": "", 
            "barcode_nested_output_folder": "0", 
            "beam_cut": "100.000000", 
            "beam_width": "32", 
            "bed_file": "", 
            "builtin_scripts": "1", 
            "calib_detect": "0", 
            "calib_max_sequence_length": "3800", 
            "calib_min_coverage": "0.600000", 
            "calib_min_sequence_length": "3000", 
            "calib_reference": "lambda_3.6kb.fasta", 
            "chunk_size": "2000", 
            "chunks_per_caller": "10000", 
            "chunks_per_runner": "208", 
            "client_id": "-1", 
            "compress_fastq": "1", 
            "cpu_threads_per_caller": "4", 
            "detect_adapter": "0", 
            "detect_barcodes": "0", 
            "detect_mid_strand_adapter": "0", 
            "detect_mid_strand_barcodes": "0", 
            "detect_primer": "0", 
            "device": "cuda:0", 
            "disable_pings": "0", 
            "disable_qscore_filtering": "0", 
            "dmean_threshold": "1.000000", 
            "dmean_win_size": "2", 
            "do_read_splitting": "0", 
            "duplex_window_size_max": "1000", 
            "duplex_window_size_min": "200", 
            "end_gap1": "40", 
            "end_gap2": "40", 
            "extend_gap1": "40", 
            "extend_gap2": "160", 
            "fast5_out": "0", 
            "flowcell": "", 
            "front_window_size": "150", 
            "gpu_runners_per_device": "12", 
            "high_priority_threshold": "10", 
            "index": "0", 
            "input_file_list": "", 
            "int8_mode": "0", 
            "jump_threshold": "1.000000", 
            "kernel_path": "", 
            "kit": "", 
            "lamp_kit": "", 
            "log_speed_frequency": "0", 
            "max_queued_reads": "2000", 
            "max_read_split_depth": "2", 
            "max_search_len": "1000", 
            "medium_priority_threshold": "4", 
            "min_length_lamp_context": "30", 
            "min_length_lamp_target": "70", 
            "min_qscore": "10.000000", 
            "min_score_adapter": "60.000000", 
            "min_score_adapter_mid": "50.000000", 
            "min_score_barcode_front": "60.000000", 
            "min_score_barcode_mask": "40.000000", 
            "min_score_barcode_mid": "50.000000", 
            "min_score_barcode_rear": "60.000000", 
            "min_score_lamp": "80.000000", 
            "min_score_lamp_mask": "50.000000", 
            "min_score_lamp_target": "50.000000", 
            "min_score_primer": "60.000000", 
            "min_score_read_splitting": "70.000000", 
            "model_file": "template_r9.4.1_450bps_sup.jsn", 
            "moves_out": "0", 
            "nested_output_folder": "0", 
            "noisiest_section_scaling_max_size": "8000", 
            "num_alignment_threads": "4", 
            "num_barcode_threads": "4", 
            "num_barcoding_buffers": "24", 
            "num_base_mod_threads": "2", 
            "num_callers": "4", 
            "num_extra_bases_trim": "0", 
            "num_mid_barcoding_buffers": "96", 
            "num_read_splitting_buffers": "16", 
            "num_read_splitting_threads": "4", 
            "num_reads_per_barcoding_buffer": "4", 
            "open_gap1": "40", 
            "open_gap2": "160", 
            "overlap": "100", 
            "override_scaling": "0", 
            "ping_segment_duration": "60", 
            "ping_url": "https://ping.oxfordnanoportal.com/basecall", 
            "post_out": "0", 
            "print_workflows": "0", 
            "progress_stats_frequency": "-1.000000", 
            "pt_median_offset": "2.500000", 
            "pt_minimum_read_start_index": "30", 
            "pt_required_adapter_drop": "30.000000", 
            "pt_scaling": "0", 
            "qscore_offset": "0.349800", 
            "qscore_scale": "0.972200", 
            "quiet": "0", 
            "read_batch_size": "4000", 
            "read_id_list": "", 
            "read_splitting_arrangement_files": "", 
            "read_splitting_score_matrix_filename": "", 
            "rear_window_size": "150", 
            "records_per_fastq": "0", 
            "recursive": "0", 
            "remora_models": "", 
            "require_barcodes_both_ends": "0", 
            "resume": "0", 
            "reverse_sequence": "0", 
            "sample_sheet": "", 
            "scaling_mad": "1.000000", 
            "scaling_med": "0.000000", 
            "score_matrix_filename": "5x5_mismatch_matrix.txt", 
            "start_gap1": "40", 
            "start_gap2": "40", 
            "stay_penalty": "1.000000", 
            "temp_bias": "1.000000", 
            "temp_weight": "1.000000", 
            "trace_categories_logs": "", 
            "trace_domains_config": "", 
            "trim_adapters": "0", 
            "trim_barcodes": "0", 
            "trim_min_events": "3", 
            "trim_primers": "0", 
            "trim_strategy": "dna", 
            "trim_threshold": "2.500000", 
            "u_substitution": "0", 
            "verbose_logs": "0"
        }, 
        "read_count": 2, 
        "reads_per_channel_dist": [
            {
                "channel": 95, 
                "count": 1
            }, 
            {
                "channel": 99, 
                "count": 1
            }
        ], 
        "run_id": "6a94eb36563a9975b38e81c41534f6590c24867c", 
        "segment_duration": 60, 
        "segment_number": 13, 
        "segment_type": "guppy-acquisition", 
        "software": {
            "analysis": "1d_basecalling", 
            "name": "guppy-basecalling", 
            "version": "6.1.2+e0556ff"
        }, 
        "tracking_id": {
            "asic_id": "3506339664", 
            "asic_id_eeprom": "5924305", 
            "asic_temp": "35.086052", 
            "asic_version": "IA02D", 
            "auto_update": "0", 
            "auto_update_source": "https://mirror.oxfordnanoportal.com/software/MinKNOW/", 
            "bream_is_standard": "0", 
            "configuration_version": "4.2.11", 
            "device_id": "MN34875", 
            "device_type": "minion", 
            "distribution_status": "stable", 
            "distribution_version": "21.02.1", 
            "exp_script_name": "sequencing_MIN106_DNA:FLO-FLG001:SQK-LSK110", 
            "exp_script_purpose": "sequencing_run", 
            "exp_start_time": "2021-08-03T19:16:19Z", 
            "flongle_adapter_id": "FA-00993", 
            "flow_cell_id": "AHC239", 
            "flow_cell_product_code": "FLO-FLG001", 
            "guppy_version": "4.3.4+ecb2805", 
            "heatsink_temp": "35.003906", 
            "hostname": "mckennalab", 
            "installation_type": "nc", 
            "local_firmware_file": "1", 
            "msg_id": "e1d2fdb6-2711-44da-82ac-d8ff71c008e7", 
            "operating_system": "ubuntu 20.04", 
            "protocol_group_id": "FEX_120", 
            "protocol_run_id": "bdcd3859-b6c1-400c-b4e6-2bc287733f75", 
            "protocol_start_time": "", 
            "protocols_version": "6.1.10", 
            "run_id": "6a94eb36563a9975b38e81c41534f6590c24867c", 
            "sample_id": "96_v2", 
            "time_stamp": "2022-06-11T02:04:21Z", 
            "usb_config": "MinION_fx3_1.1.1_ONT#MinION_fpga_1.1.0#bulk#Auto", 
            "version": "4.2.5"
        }
    }, 
    {
        "aggregation": "segment", 
        "analysis_id": "e2692740-f853-4dfa-99c2-3309f56cdfdb", 
        "basecall_1d": {
            "exit_status_dist": {
                "pass": 4
            }, 
            "qscore_dist_temp": [
                {
                    "count": 2, 
                    "mean_qscore": 13.0
                }, 
                {
                    "count": 2, 
                    "mean_qscore": 14.0
                }
            ], 
            "qscore_sum_temp": {
                "count": 4, 
                "mean": 13.804849624633789, 
                "sum": 55.219398498535156
            }, 
            "read_len_events_sum_temp": 106141, 
            "seq_len_bases_dist_temp": [
                {
                    "count": 4, 
                    "length": 0.0
                }
            ], 
            "seq_len_bases_sum_temp": 4, 
            "seq_len_events_dist_temp": [
                {
                    "count": 1, 
                    "length": 18000.0
                }, 
                {
                    "count": 1, 
                    "length": 23000.0
                }, 
                {
                    "count": 1, 
                    "length": 29000.0
                }, 
                {
                    "count": 1, 
                    "length": 34000.0
                }
            ], 
            "speed_bases_per_second_dist_temp": [
                {
                    "count": 4, 
                    "speed": 1.0
                }
            ], 
            "strand_median_pa": {
                "count": 4, 
                "mean": 69.14209747314453, 
                "sum": 276.5683898925781
            }, 
            "strand_sd_pa": {
                "count": 4, 
                "mean": 9.549036979675293, 
                "sum": 38.19614791870117
            }
        }, 
        "channel_count": 4, 
        "context_tags": {
            "barcoding_enabled": "0", 
            "experiment_duration_set": "1680", 
            "experiment_type": "genomic_dna", 
            "local_basecalling": "0", 
            "package": "bream4", 
            "package_version": "6.1.10", 
            "sample_frequency": "4000", 
            "sequencing_kit": "sqk-lsk110"
        }, 
        "latest_run_time": 50318.1796875, 
        "levels_sums": {
            "count": 4, 
            "mean": 188.97772216796875, 
            "open_pore_level_sum": 755.910888671875
        }, 
        "opts": {
            "adapter_pt_range_scale": "5.200000", 
            "additional_lamp_context_bases": "2", 
            "align_ref": "", 
            "align_type": "auto", 
            "allow_inferior_barcodes": "0", 
            "as_cpu_threads_per_scaler": "2", 
            "as_gpu_runners_per_device": "2", 
            "as_model_file": "", 
            "as_num_scalers": "4", 
            "as_reads_per_runner": "32", 
            "bam_methylation_threshold": "5.000000", 
            "bam_out": "0", 
            "barcode_kits": "", 
            "barcode_nested_output_folder": "0", 
            "beam_cut": "100.000000", 
            "beam_width": "32", 
            "bed_file": "", 
            "builtin_scripts": "1", 
            "calib_detect": "0", 
            "calib_max_sequence_length": "3800", 
            "calib_min_coverage": "0.600000", 
            "calib_min_sequence_length": "3000", 
            "calib_reference": "lambda_3.6kb.fasta", 
            "chunk_size": "2000", 
            "chunks_per_caller": "10000", 
            "chunks_per_runner": "208", 
            "client_id": "-1", 
            "compress_fastq": "1", 
            "cpu_threads_per_caller": "4", 
            "detect_adapter": "0", 
            "detect_barcodes": "0", 
            "detect_mid_strand_adapter": "0", 
            "detect_mid_strand_barcodes": "0", 
            "detect_primer": "0", 
            "device": "cuda:0", 
            "disable_pings": "0", 
            "disable_qscore_filtering": "0", 
            "dmean_threshold": "1.000000", 
            "dmean_win_size": "2", 
            "do_read_splitting": "0", 
            "duplex_window_size_max": "1000", 
            "duplex_window_size_min": "200", 
            "end_gap1": "40", 
            "end_gap2": "40", 
            "extend_gap1": "40", 
            "extend_gap2": "160", 
            "fast5_out": "0", 
            "flowcell": "", 
            "front_window_size": "150", 
            "gpu_runners_per_device": "12", 
            "high_priority_threshold": "10", 
            "index": "0", 
            "input_file_list": "", 
            "int8_mode": "0", 
            "jump_threshold": "1.000000", 
            "kernel_path": "", 
            "kit": "", 
            "lamp_kit": "", 
            "log_speed_frequency": "0", 
            "max_queued_reads": "2000", 
            "max_read_split_depth": "2", 
            "max_search_len": "1000", 
            "medium_priority_threshold": "4", 
            "min_length_lamp_context": "30", 
            "min_length_lamp_target": "70", 
            "min_qscore": "10.000000", 
            "min_score_adapter": "60.000000", 
            "min_score_adapter_mid": "50.000000", 
            "min_score_barcode_front": "60.000000", 
            "min_score_barcode_mask": "40.000000", 
            "min_score_barcode_mid": "50.000000", 
            "min_score_barcode_rear": "60.000000", 
            "min_score_lamp": "80.000000", 
            "min_score_lamp_mask": "50.000000", 
            "min_score_lamp_target": "50.000000", 
            "min_score_primer": "60.000000", 
            "min_score_read_splitting": "70.000000", 
            "model_file": "template_r9.4.1_450bps_sup.jsn", 
            "moves_out": "0", 
            "nested_output_folder": "0", 
            "noisiest_section_scaling_max_size": "8000", 
            "num_alignment_threads": "4", 
            "num_barcode_threads": "4", 
            "num_barcoding_buffers": "24", 
            "num_base_mod_threads": "2", 
            "num_callers": "4", 
            "num_extra_bases_trim": "0", 
            "num_mid_barcoding_buffers": "96", 
            "num_read_splitting_buffers": "16", 
            "num_read_splitting_threads": "4", 
            "num_reads_per_barcoding_buffer": "4", 
            "open_gap1": "40", 
            "open_gap2": "160", 
            "overlap": "100", 
            "override_scaling": "0", 
            "ping_segment_duration": "60", 
            "ping_url": "https://ping.oxfordnanoportal.com/basecall", 
            "post_out": "0", 
            "print_workflows": "0", 
            "progress_stats_frequency": "-1.000000", 
            "pt_median_offset": "2.500000", 
            "pt_minimum_read_start_index": "30", 
            "pt_required_adapter_drop": "30.000000", 
            "pt_scaling": "0", 
            "qscore_offset": "0.349800", 
            "qscore_scale": "0.972200", 
            "quiet": "0", 
            "read_batch_size": "4000", 
            "read_id_list": "", 
            "read_splitting_arrangement_files": "", 
            "read_splitting_score_matrix_filename": "", 
            "rear_window_size": "150", 
            "records_per_fastq": "0", 
            "recursive": "0", 
            "remora_models": "", 
            "require_barcodes_both_ends": "0", 
            "resume": "0", 
            "reverse_sequence": "0", 
            "sample_sheet": "", 
            "scaling_mad": "1.000000", 
            "scaling_med": "0.000000", 
            "score_matrix_filename": "5x5_mismatch_matrix.txt", 
            "start_gap1": "40", 
            "start_gap2": "40", 
            "stay_penalty": "1.000000", 
            "temp_bias": "1.000000", 
            "temp_weight": "1.000000", 
            "trace_categories_logs": "", 
            "trace_domains_config": "", 
            "trim_adapters": "0", 
            "trim_barcodes": "0", 
            "trim_min_events": "3", 
            "trim_primers": "0", 
            "trim_strategy": "dna", 
            "trim_threshold": "2.500000", 
            "u_substitution": "0", 
            "verbose_logs": "0"
        }, 
        "read_count": 4, 
        "reads_per_channel_dist": [
            {
                "channel": 29, 
                "count": 1
            }, 
            {
                "channel": 36, 
                "count": 1
            }, 
            {
                "channel": 71, 
                "count": 1
            }, 
            {
                "channel": 76, 
                "count": 1
            }
        ], 
        "run_id": "6a94eb36563a9975b38e81c41534f6590c24867c", 
        "segment_duration": 60, 
        "segment_number": 14, 
        "segment_type": "guppy-acquisition", 
        "software": {
            "analysis": "1d_basecalling", 
            "name": "guppy-basecalling", 
            "version": "6.1.2+e0556ff"
        }, 
        "tracking_id": {
            "asic_id": "3506339664", 
            "asic_id_eeprom": "5924305", 
            "asic_temp": "35.086052", 
            "asic_version": "IA02D", 
            "auto_update": "0", 
            "auto_update_source": "https://mirror.oxfordnanoportal.com/software/MinKNOW/", 
            "bream_is_standard": "0", 
            "configuration_version": "4.2.11", 
            "device_id": "MN34875", 
            "device_type": "minion", 
            "distribution_status": "stable", 
            "distribution_version": "21.02.1", 
            "exp_script_name": "sequencing_MIN106_DNA:FLO-FLG001:SQK-LSK110", 
            "exp_script_purpose": "sequencing_run", 
            "exp_start_time": "2021-08-03T19:16:19Z", 
            "flongle_adapter_id": "FA-00993", 
            "flow_cell_id": "AHC239", 
            "flow_cell_product_code": "FLO-FLG001", 
            "guppy_version": "4.3.4+ecb2805", 
            "heatsink_temp": "35.003906", 
            "hostname": "mckennalab", 
            "installation_type": "nc", 
            "local_firmware_file": "1", 
            "msg_id": "de153322-6c4a-4aca-9eb7-9499830f3465", 
            "operating_system": "ubuntu 20.04", 
            "protocol_group_id": "FEX_120", 
            "protocol_run_id": "bdcd3859-b6c1-400c-b4e6-2bc287733f75", 
            "protocol_start_time": "", 
            "protocols_version": "6.1.10", 
            "run_id": "6a94eb36563a9975b38e81c41534f6590c24867c", 
            "sample_id": "96_v2", 
            "time_stamp": "2022-06-11T02:04:21Z", 
            "usb_config": "MinION_fx3_1.1.1_ONT#MinION_fpga_1.1.0#bulk#Auto", 
            "version": "4.2.5"
        }
    }, 
    {
        "aggregation": "segment", 
        "analysis_id": "e2692740-f853-4dfa-99c2-3309f56cdfdb", 
        "basecall_1d": {
            "exit_status_dist": {
                "pass": 3
            }, 
            "qscore_dist_temp": [
                {
                    "count": 1, 
                    "mean_qscore": 11.5
                }, 
                {
                    "count": 1, 
                    "mean_qscore": 14.0
                }, 
                {
                    "count": 1, 
                    "mean_qscore": 14.5
                }
            ], 
            "qscore_sum_temp": {
                "count": 3, 
                "mean": 13.500847816467285, 
                "sum": 40.50254440307617
            }, 
            "read_len_events_sum_temp": 90735, 
            "seq_len_bases_dist_temp": [
                {
                    "count": 3, 
                    "length": 0.0
                }
            ], 
            "seq_len_bases_sum_temp": 3, 
            "seq_len_events_dist_temp": [
                {
                    "count": 1, 
                    "length": 24000.0
                }, 
                {
                    "count": 1, 
                    "length": 30000.0
                }, 
                {
                    "count": 1, 
                    "length": 35000.0
                }
            ], 
            "speed_bases_per_second_dist_temp": [
                {
                    "count": 3, 
                    "speed": 1.0
                }
            ], 
            "strand_median_pa": {
                "count": 3, 
                "mean": 72.26616668701172, 
                "sum": 216.79849243164062
            }, 
            "strand_sd_pa": {
                "count": 3, 
                "mean": 9.725871086120605, 
                "sum": 29.1776123046875
            }
        }, 
        "channel_count": 3, 
        "context_tags": {
            "barcoding_enabled": "0", 
            "experiment_duration_set": "1680", 
            "experiment_type": "genomic_dna", 
            "local_basecalling": "0", 
            "package": "bream4", 
            "package_version": "6.1.10", 
            "sample_frequency": "4000", 
            "sequencing_kit": "sqk-lsk110"
        }, 
        "latest_run_time": 58619.8359375, 
        "levels_sums": {
            "count": 3, 
            "mean": 204.9744415283203, 
            "open_pore_level_sum": 614.92333984375
        }, 
        "opts": {
            "adapter_pt_range_scale": "5.200000", 
            "additional_lamp_context_bases": "2", 
            "align_ref": "", 
            "align_type": "auto", 
            "allow_inferior_barcodes": "0", 
            "as_cpu_threads_per_scaler": "2", 
            "as_gpu_runners_per_device": "2", 
            "as_model_file": "", 
            "as_num_scalers": "4", 
            "as_reads_per_runner": "32", 
            "bam_methylation_threshold": "5.000000", 
            "bam_out": "0", 
            "barcode_kits": "", 
            "barcode_nested_output_folder": "0", 
            "beam_cut": "100.000000", 
            "beam_width": "32", 
            "bed_file": "", 
            "builtin_scripts": "1", 
            "calib_detect": "0", 
            "calib_max_sequence_length": "3800", 
            "calib_min_coverage": "0.600000", 
            "calib_min_sequence_length": "3000", 
            "calib_reference": "lambda_3.6kb.fasta", 
            "chunk_size": "2000", 
            "chunks_per_caller": "10000", 
            "chunks_per_runner": "208", 
            "client_id": "-1", 
            "compress_fastq": "1", 
            "cpu_threads_per_caller": "4", 
            "detect_adapter": "0", 
            "detect_barcodes": "0", 
            "detect_mid_strand_adapter": "0", 
            "detect_mid_strand_barcodes": "0", 
            "detect_primer": "0", 
            "device": "cuda:0", 
            "disable_pings": "0", 
            "disable_qscore_filtering": "0", 
            "dmean_threshold": "1.000000", 
            "dmean_win_size": "2", 
            "do_read_splitting": "0", 
            "duplex_window_size_max": "1000", 
            "duplex_window_size_min": "200", 
            "end_gap1": "40", 
            "end_gap2": "40", 
            "extend_gap1": "40", 
            "extend_gap2": "160", 
            "fast5_out": "0", 
            "flowcell": "", 
            "front_window_size": "150", 
            "gpu_runners_per_device": "12", 
            "high_priority_threshold": "10", 
            "index": "0", 
            "input_file_list": "", 
            "int8_mode": "0", 
            "jump_threshold": "1.000000", 
            "kernel_path": "", 
            "kit": "", 
            "lamp_kit": "", 
            "log_speed_frequency": "0", 
            "max_queued_reads": "2000", 
            "max_read_split_depth": "2", 
            "max_search_len": "1000", 
            "medium_priority_threshold": "4", 
            "min_length_lamp_context": "30", 
            "min_length_lamp_target": "70", 
            "min_qscore": "10.000000", 
            "min_score_adapter": "60.000000", 
            "min_score_adapter_mid": "50.000000", 
            "min_score_barcode_front": "60.000000", 
            "min_score_barcode_mask": "40.000000", 
            "min_score_barcode_mid": "50.000000", 
            "min_score_barcode_rear": "60.000000", 
            "min_score_lamp": "80.000000", 
            "min_score_lamp_mask": "50.000000", 
            "min_score_lamp_target": "50.000000", 
            "min_score_primer": "60.000000", 
            "min_score_read_splitting": "70.000000", 
            "model_file": "template_r9.4.1_450bps_sup.jsn", 
            "moves_out": "0", 
            "nested_output_folder": "0", 
            "noisiest_section_scaling_max_size": "8000", 
            "num_alignment_threads": "4", 
            "num_barcode_threads": "4", 
            "num_barcoding_buffers": "24", 
            "num_base_mod_threads": "2", 
            "num_callers": "4", 
            "num_extra_bases_trim": "0", 
            "num_mid_barcoding_buffers": "96", 
            "num_read_splitting_buffers": "16", 
            "num_read_splitting_threads": "4", 
            "num_reads_per_barcoding_buffer": "4", 
            "open_gap1": "40", 
            "open_gap2": "160", 
            "overlap": "100", 
            "override_scaling": "0", 
            "ping_segment_duration": "60", 
            "ping_url": "https://ping.oxfordnanoportal.com/basecall", 
            "post_out": "0", 
            "print_workflows": "0", 
            "progress_stats_frequency": "-1.000000", 
            "pt_median_offset": "2.500000", 
            "pt_minimum_read_start_index": "30", 
            "pt_required_adapter_drop": "30.000000", 
            "pt_scaling": "0", 
            "qscore_offset": "0.349800", 
            "qscore_scale": "0.972200", 
            "quiet": "0", 
            "read_batch_size": "4000", 
            "read_id_list": "", 
            "read_splitting_arrangement_files": "", 
            "read_splitting_score_matrix_filename": "", 
            "rear_window_size": "150", 
            "records_per_fastq": "0", 
            "recursive": "0", 
            "remora_models": "", 
            "require_barcodes_both_ends": "0", 
            "resume": "0", 
            "reverse_sequence": "0", 
            "sample_sheet": "", 
            "scaling_mad": "1.000000", 
            "scaling_med": "0.000000", 
            "score_matrix_filename": "5x5_mismatch_matrix.txt", 
            "start_gap1": "40", 
            "start_gap2": "40", 
            "stay_penalty": "1.000000", 
            "temp_bias": "1.000000", 
            "temp_weight": "1.000000", 
            "trace_categories_logs": "", 
            "trace_domains_config": "", 
            "trim_adapters": "0", 
            "trim_barcodes": "0", 
            "trim_min_events": "3", 
            "trim_primers": "0", 
            "trim_strategy": "dna", 
            "trim_threshold": "2.500000", 
            "u_substitution": "0", 
            "verbose_logs": "0"
        }, 
        "read_count": 3, 
        "reads_per_channel_dist": [
            {
                "channel": 16, 
                "count": 1
            }, 
            {
                "channel": 31, 
                "count": 1
            }, 
            {
                "channel": 35, 
                "count": 1
            }
        ], 
        "run_id": "6a94eb36563a9975b38e81c41534f6590c24867c", 
        "segment_duration": 60, 
        "segment_number": 17, 
        "segment_type": "guppy-acquisition", 
        "software": {
            "analysis": "1d_basecalling", 
            "name": "guppy-basecalling", 
            "version": "6.1.2+e0556ff"
        }, 
        "tracking_id": {
            "asic_id": "3506339664", 
            "asic_id_eeprom": "5924305", 
            "asic_temp": "35.086052", 
            "asic_version": "IA02D", 
            "auto_update": "0", 
            "auto_update_source": "https://mirror.oxfordnanoportal.com/software/MinKNOW/", 
            "bream_is_standard": "0", 
            "configuration_version": "4.2.11", 
            "device_id": "MN34875", 
            "device_type": "minion", 
            "distribution_status": "stable", 
            "distribution_version": "21.02.1", 
            "exp_script_name": "sequencing_MIN106_DNA:FLO-FLG001:SQK-LSK110", 
            "exp_script_purpose": "sequencing_run", 
            "exp_start_time": "2021-08-03T19:16:19Z", 
            "flongle_adapter_id": "FA-00993", 
            "flow_cell_id": "AHC239", 
            "flow_cell_product_code": "FLO-FLG001", 
            "guppy_version": "4.3.4+ecb2805", 
            "heatsink_temp": "35.003906", 
            "hostname": "mckennalab", 
            "installation_type": "nc", 
            "local_firmware_file": "1", 
            "msg_id": "93553e33-c970-4126-aa46-4f7371bb06f4", 
            "operating_system": "ubuntu 20.04", 
            "protocol_group_id": "FEX_120", 
            "protocol_run_id": "bdcd3859-b6c1-400c-b4e6-2bc287733f75", 
            "protocol_start_time": "", 
            "protocols_version": "6.1.10", 
            "run_id": "6a94eb36563a9975b38e81c41534f6590c24867c", 
            "sample_id": "96_v2", 
            "time_stamp": "2022-06-11T02:04:21Z", 
            "usb_config": "MinION_fx3_1.1.1_ONT#MinION_fpga_1.1.0#bulk#Auto", 
            "version": "4.2.5"
        }
    }, 
    {
        "aggregation": "segment", 
        "analysis_id": "e2692740-f853-4dfa-99c2-3309f56cdfdb", 
        "basecall_1d": {
            "exit_status_dist": {
                "pass": 1
            }, 
            "qscore_dist_temp": [
                {
                    "count": 1, 
                    "mean_qscore": 13.5
                }
            ], 
            "qscore_sum_temp": {
                "count": 1, 
                "mean": 13.771560668945312, 
                "sum": 13.771560668945312
            }, 
            "read_len_events_sum_temp": 22482, 
            "seq_len_bases_dist_temp": [
                {
                    "count": 1, 
                    "length": 0.0
                }
            ], 
            "seq_len_bases_sum_temp": 1, 
            "seq_len_events_dist_temp": [
                {
                    "count": 1, 
                    "length": 22000.0
                }
            ], 
            "speed_bases_per_second_dist_temp": [
                {
                    "count": 1, 
                    "speed": 1.0
                }
            ], 
            "strand_median_pa": {
                "count": 1, 
                "mean": 68.96526336669922, 
                "sum": 68.96526336669922
            }, 
            "strand_sd_pa": {
                "count": 1, 
                "mean": 9.37220287322998, 
                "sum": 9.37220287322998
            }
        }, 
        "channel_count": 1, 
        "context_tags": {
            "barcoding_enabled": "0", 
            "experiment_duration_set": "1680", 
            "experiment_type": "genomic_dna", 
            "local_basecalling": "0", 
            "package": "bream4", 
            "package_version": "6.1.10", 
            "sample_frequency": "4000", 
            "sequencing_kit": "sqk-lsk110"
        }, 
        "latest_run_time": 62672.2734375, 
        "levels_sums": {
            "count": 1, 
            "mean": 69.42503356933594, 
            "open_pore_level_sum": 69.42503356933594
        }, 
        "opts": {
            "adapter_pt_range_scale": "5.200000", 
            "additional_lamp_context_bases": "2", 
            "align_ref": "", 
            "align_type": "auto", 
            "allow_inferior_barcodes": "0", 
            "as_cpu_threads_per_scaler": "2", 
            "as_gpu_runners_per_device": "2", 
            "as_model_file": "", 
            "as_num_scalers": "4", 
            "as_reads_per_runner": "32", 
            "bam_methylation_threshold": "5.000000", 
            "bam_out": "0", 
            "barcode_kits": "", 
            "barcode_nested_output_folder": "0", 
            "beam_cut": "100.000000", 
            "beam_width": "32", 
            "bed_file": "", 
            "builtin_scripts": "1", 
            "calib_detect": "0", 
            "calib_max_sequence_length": "3800", 
            "calib_min_coverage": "0.600000", 
            "calib_min_sequence_length": "3000", 
            "calib_reference": "lambda_3.6kb.fasta", 
            "chunk_size": "2000", 
            "chunks_per_caller": "10000", 
            "chunks_per_runner": "208", 
            "client_id": "-1", 
            "compress_fastq": "1", 
            "cpu_threads_per_caller": "4", 
            "detect_adapter": "0", 
            "detect_barcodes": "0", 
            "detect_mid_strand_adapter": "0", 
            "detect_mid_strand_barcodes": "0", 
            "detect_primer": "0", 
            "device": "cuda:0", 
            "disable_pings": "0", 
            "disable_qscore_filtering": "0", 
            "dmean_threshold": "1.000000", 
            "dmean_win_size": "2", 
            "do_read_splitting": "0", 
            "duplex_window_size_max": "1000", 
            "duplex_window_size_min": "200", 
            "end_gap1": "40", 
            "end_gap2": "40", 
            "extend_gap1": "40", 
            "extend_gap2": "160", 
            "fast5_out": "0", 
            "flowcell": "", 
            "front_window_size": "150", 
            "gpu_runners_per_device": "12", 
            "high_priority_threshold": "10", 
            "index": "0", 
            "input_file_list": "", 
            "int8_mode": "0", 
            "jump_threshold": "1.000000", 
            "kernel_path": "", 
            "kit": "", 
            "lamp_kit": "", 
            "log_speed_frequency": "0", 
            "max_queued_reads": "2000", 
            "max_read_split_depth": "2", 
            "max_search_len": "1000", 
            "medium_priority_threshold": "4", 
            "min_length_lamp_context": "30", 
            "min_length_lamp_target": "70", 
            "min_qscore": "10.000000", 
            "min_score_adapter": "60.000000", 
            "min_score_adapter_mid": "50.000000", 
            "min_score_barcode_front": "60.000000", 
            "min_score_barcode_mask": "40.000000", 
            "min_score_barcode_mid": "50.000000", 
            "min_score_barcode_rear": "60.000000", 
            "min_score_lamp": "80.000000", 
            "min_score_lamp_mask": "50.000000", 
            "min_score_lamp_target": "50.000000", 
            "min_score_primer": "60.000000", 
            "min_score_read_splitting": "70.000000", 
            "model_file": "template_r9.4.1_450bps_sup.jsn", 
            "moves_out": "0", 
            "nested_output_folder": "0", 
            "noisiest_section_scaling_max_size": "8000", 
            "num_alignment_threads": "4", 
            "num_barcode_threads": "4", 
            "num_barcoding_buffers": "24", 
            "num_base_mod_threads": "2", 
            "num_callers": "4", 
            "num_extra_bases_trim": "0", 
            "num_mid_barcoding_buffers": "96", 
            "num_read_splitting_buffers": "16", 
            "num_read_splitting_threads": "4", 
            "num_reads_per_barcoding_buffer": "4", 
            "open_gap1": "40", 
            "open_gap2": "160", 
            "overlap": "100", 
            "override_scaling": "0", 
            "ping_segment_duration": "60", 
            "ping_url": "https://ping.oxfordnanoportal.com/basecall", 
            "post_out": "0", 
            "print_workflows": "0", 
            "progress_stats_frequency": "-1.000000", 
            "pt_median_offset": "2.500000", 
            "pt_minimum_read_start_index": "30", 
            "pt_required_adapter_drop": "30.000000", 
            "pt_scaling": "0", 
            "qscore_offset": "0.349800", 
            "qscore_scale": "0.972200", 
            "quiet": "0", 
            "read_batch_size": "4000", 
            "read_id_list": "", 
            "read_splitting_arrangement_files": "", 
            "read_splitting_score_matrix_filename": "", 
            "rear_window_size": "150", 
            "records_per_fastq": "0", 
            "recursive": "0", 
            "remora_models": "", 
            "require_barcodes_both_ends": "0", 
            "resume": "0", 
            "reverse_sequence": "0", 
            "sample_sheet": "", 
            "scaling_mad": "1.000000", 
            "scaling_med": "0.000000", 
            "score_matrix_filename": "5x5_mismatch_matrix.txt", 
            "start_gap1": "40", 
            "start_gap2": "40", 
            "stay_penalty": "1.000000", 
            "temp_bias": "1.000000", 
            "temp_weight": "1.000000", 
            "trace_categories_logs": "", 
            "trace_domains_config": "", 
            "trim_adapters": "0", 
            "trim_barcodes": "0", 
            "trim_min_events": "3", 
            "trim_primers": "0", 
            "trim_strategy": "dna", 
            "trim_threshold": "2.500000", 
            "u_substitution": "0", 
            "verbose_logs": "0"
        }, 
        "read_count": 1, 
        "reads_per_channel_dist": [
            {
                "channel": 111, 
                "count": 1
            }
        ], 
        "run_id": "6a94eb36563a9975b38e81c41534f6590c24867c", 
        "segment_duration": 60, 
        "segment_number": 18, 
        "segment_type": "guppy-acquisition", 
        "software": {
            "analysis": "1d_basecalling", 
            "name": "guppy-basecalling", 
            "version": "6.1.2+e0556ff"
        }, 
        "tracking_id": {
            "asic_id": "3506339664", 
            "asic_id_eeprom": "5924305", 
            "asic_temp": "35.086052", 
            "asic_version": "IA02D", 
            "auto_update": "0", 
            "auto_update_source": "https://mirror.oxfordnanoportal.com/software/MinKNOW/", 
            "bream_is_standard": "0", 
            "configuration_version": "4.2.11", 
            "device_id": "MN34875", 
            "device_type": "minion", 
            "distribution_status": "stable", 
            "distribution_version": "21.02.1", 
            "exp_script_name": "sequencing_MIN106_DNA:FLO-FLG001:SQK-LSK110", 
            "exp_script_purpose": "sequencing_run", 
            "exp_start_time": "2021-08-03T19:16:19Z", 
            "flongle_adapter_id": "FA-00993", 
            "flow_cell_id": "AHC239", 
            "flow_cell_product_code": "FLO-FLG001", 
            "guppy_version": "4.3.4+ecb2805", 
            "heatsink_temp": "35.003906", 
            "hostname": "mckennalab", 
            "installation_type": "nc", 
            "local_firmware_file": "1", 
            "msg_id": "ad79c373-aeda-4a8d-bf21-0d48ba5a5d44", 
            "operating_system": "ubuntu 20.04", 
            "protocol_group_id": "FEX_120", 
            "protocol_run_id": "bdcd3859-b6c1-400c-b4e6-2bc287733f75", 
            "protocol_start_time": "", 
            "protocols_version": "6.1.10", 
            "run_id": "6a94eb36563a9975b38e81c41534f6590c24867c", 
            "sample_id": "96_v2", 
            "time_stamp": "2022-06-11T02:04:21Z", 
            "usb_config": "MinION_fx3_1.1.1_ONT#MinION_fpga_1.1.0#bulk#Auto", 
            "version": "4.2.5"
        }
    }, 
    {
        "aggregation": "segment", 
        "analysis_id": "e2692740-f853-4dfa-99c2-3309f56cdfdb", 
        "basecall_1d": {
            "exit_status_dist": {
                "pass": 1
            }, 
            "qscore_dist_temp": [
                {
                    "count": 1, 
                    "mean_qscore": 13.5
                }
            ], 
            "qscore_sum_temp": {
                "count": 1, 
                "mean": 13.938974380493164, 
                "sum": 13.938974380493164
            }, 
            "read_len_events_sum_temp": 33056, 
            "seq_len_bases_dist_temp": [
                {
                    "count": 1, 
                    "length": 0.0
                }
            ], 
            "seq_len_bases_sum_temp": 1, 
            "seq_len_events_dist_temp": [
                {
                    "count": 1, 
                    "length": 33000.0
                }
            ], 
            "speed_bases_per_second_dist_temp": [
                {
                    "count": 1, 
                    "speed": 1.0
                }
            ], 
            "strand_median_pa": {
                "count": 1, 
                "mean": 69.67259979248047, 
                "sum": 69.67259979248047
            }, 
            "strand_sd_pa": {
                "count": 1, 
                "mean": 9.018534660339355, 
                "sum": 9.018534660339355
            }
        }, 
        "channel_count": 1, 
        "context_tags": {
            "barcoding_enabled": "0", 
            "experiment_duration_set": "1680", 
            "experiment_type": "genomic_dna", 
            "local_basecalling": "0", 
            "package": "bream4", 
            "package_version": "6.1.10", 
            "sample_frequency": "4000", 
            "sequencing_kit": "sqk-lsk110"
        }, 
        "latest_run_time": 65498.953125, 
        "levels_sums": {
            "count": 1, 
            "mean": 198.4498748779297, 
            "open_pore_level_sum": 198.4498748779297
        }, 
        "opts": {
            "adapter_pt_range_scale": "5.200000", 
            "additional_lamp_context_bases": "2", 
            "align_ref": "", 
            "align_type": "auto", 
            "allow_inferior_barcodes": "0", 
            "as_cpu_threads_per_scaler": "2", 
            "as_gpu_runners_per_device": "2", 
            "as_model_file": "", 
            "as_num_scalers": "4", 
            "as_reads_per_runner": "32", 
            "bam_methylation_threshold": "5.000000", 
            "bam_out": "0", 
            "barcode_kits": "", 
            "barcode_nested_output_folder": "0", 
            "beam_cut": "100.000000", 
            "beam_width": "32", 
            "bed_file": "", 
            "builtin_scripts": "1", 
            "calib_detect": "0", 
            "calib_max_sequence_length": "3800", 
            "calib_min_coverage": "0.600000", 
            "calib_min_sequence_length": "3000", 
            "calib_reference": "lambda_3.6kb.fasta", 
            "chunk_size": "2000", 
            "chunks_per_caller": "10000", 
            "chunks_per_runner": "208", 
            "client_id": "-1", 
            "compress_fastq": "1", 
            "cpu_threads_per_caller": "4", 
            "detect_adapter": "0", 
            "detect_barcodes": "0", 
            "detect_mid_strand_adapter": "0", 
            "detect_mid_strand_barcodes": "0", 
            "detect_primer": "0", 
            "device": "cuda:0", 
            "disable_pings": "0", 
            "disable_qscore_filtering": "0", 
            "dmean_threshold": "1.000000", 
            "dmean_win_size": "2", 
            "do_read_splitting": "0", 
            "duplex_window_size_max": "1000", 
            "duplex_window_size_min": "200", 
            "end_gap1": "40", 
            "end_gap2": "40", 
            "extend_gap1": "40", 
            "extend_gap2": "160", 
            "fast5_out": "0", 
            "flowcell": "", 
            "front_window_size": "150", 
            "gpu_runners_per_device": "12", 
            "high_priority_threshold": "10", 
            "index": "0", 
            "input_file_list": "", 
            "int8_mode": "0", 
            "jump_threshold": "1.000000", 
            "kernel_path": "", 
            "kit": "", 
            "lamp_kit": "", 
            "log_speed_frequency": "0", 
            "max_queued_reads": "2000", 
            "max_read_split_depth": "2", 
            "max_search_len": "1000", 
            "medium_priority_threshold": "4", 
            "min_length_lamp_context": "30", 
            "min_length_lamp_target": "70", 
            "min_qscore": "10.000000", 
            "min_score_adapter": "60.000000", 
            "min_score_adapter_mid": "50.000000", 
            "min_score_barcode_front": "60.000000", 
            "min_score_barcode_mask": "40.000000", 
            "min_score_barcode_mid": "50.000000", 
            "min_score_barcode_rear": "60.000000", 
            "min_score_lamp": "80.000000", 
            "min_score_lamp_mask": "50.000000", 
            "min_score_lamp_target": "50.000000", 
            "min_score_primer": "60.000000", 
            "min_score_read_splitting": "70.000000", 
            "model_file": "template_r9.4.1_450bps_sup.jsn", 
            "moves_out": "0", 
            "nested_output_folder": "0", 
            "noisiest_section_scaling_max_size": "8000", 
            "num_alignment_threads": "4", 
            "num_barcode_threads": "4", 
            "num_barcoding_buffers": "24", 
            "num_base_mod_threads": "2", 
            "num_callers": "4", 
            "num_extra_bases_trim": "0", 
            "num_mid_barcoding_buffers": "96", 
            "num_read_splitting_buffers": "16", 
            "num_read_splitting_threads": "4", 
            "num_reads_per_barcoding_buffer": "4", 
            "open_gap1": "40", 
            "open_gap2": "160", 
            "overlap": "100", 
            "override_scaling": "0", 
            "ping_segment_duration": "60", 
            "ping_url": "https://ping.oxfordnanoportal.com/basecall", 
            "post_out": "0", 
            "print_workflows": "0", 
            "progress_stats_frequency": "-1.000000", 
            "pt_median_offset": "2.500000", 
            "pt_minimum_read_start_index": "30", 
            "pt_required_adapter_drop": "30.000000", 
            "pt_scaling": "0", 
            "qscore_offset": "0.349800", 
            "qscore_scale": "0.972200", 
            "quiet": "0", 
            "read_batch_size": "4000", 
            "read_id_list": "", 
            "read_splitting_arrangement_files": "", 
            "read_splitting_score_matrix_filename": "", 
            "rear_window_size": "150", 
            "records_per_fastq": "0", 
            "recursive": "0", 
            "remora_models": "", 
            "require_barcodes_both_ends": "0", 
            "resume": "0", 
            "reverse_sequence": "0", 
            "sample_sheet": "", 
            "scaling_mad": "1.000000", 
            "scaling_med": "0.000000", 
            "score_matrix_filename": "5x5_mismatch_matrix.txt", 
            "start_gap1": "40", 
            "start_gap2": "40", 
            "stay_penalty": "1.000000", 
            "temp_bias": "1.000000", 
            "temp_weight": "1.000000", 
            "trace_categories_logs": "", 
            "trace_domains_config": "", 
            "trim_adapters": "0", 
            "trim_barcodes": "0", 
            "trim_min_events": "3", 
            "trim_primers": "0", 
            "trim_strategy": "dna", 
            "trim_threshold": "2.500000", 
            "u_substitution": "0", 
            "verbose_logs": "0"
        }, 
        "read_count": 1, 
        "reads_per_channel_dist": [
            {
                "channel": 29, 
                "count": 1
            }
        ], 
        "run_id": "6a94eb36563a9975b38e81c41534f6590c24867c", 
        "segment_duration": 60, 
        "segment_number": 19, 
        "segment_type": "guppy-acquisition", 
        "software": {
            "analysis": "1d_basecalling", 
            "name": "guppy-basecalling", 
            "version": "6.1.2+e0556ff"
        }, 
        "tracking_id": {
            "asic_id": "3506339664", 
            "asic_id_eeprom": "5924305", 
            "asic_temp": "35.086052", 
            "asic_version": "IA02D", 
            "auto_update": "0", 
            "auto_update_source": "https://mirror.oxfordnanoportal.com/software/MinKNOW/", 
            "bream_is_standard": "0", 
            "configuration_version": "4.2.11", 
            "device_id": "MN34875", 
            "device_type": "minion", 
            "distribution_status": "stable", 
            "distribution_version": "21.02.1", 
            "exp_script_name": "sequencing_MIN106_DNA:FLO-FLG001:SQK-LSK110", 
            "exp_script_purpose": "sequencing_run", 
            "exp_start_time": "2021-08-03T19:16:19Z", 
            "flongle_adapter_id": "FA-00993", 
            "flow_cell_id": "AHC239", 
            "flow_cell_product_code": "FLO-FLG001", 
            "guppy_version": "4.3.4+ecb2805", 
            "heatsink_temp": "35.003906", 
            "hostname": "mckennalab", 
            "installation_type": "nc", 
            "local_firmware_file": "1", 
            "msg_id": "7682dfd5-fad8-4ce6-8808-767a6384f229", 
            "operating_system": "ubuntu 20.04", 
            "protocol_group_id": "FEX_120", 
            "protocol_run_id": "bdcd3859-b6c1-400c-b4e6-2bc287733f75", 
            "protocol_start_time": "", 
            "protocols_version": "6.1.10", 
            "run_id": "6a94eb36563a9975b38e81c41534f6590c24867c", 
            "sample_id": "96_v2", 
            "time_stamp": "2022-06-11T02:04:21Z", 
            "usb_config": "MinION_fx3_1.1.1_ONT#MinION_fpga_1.1.0#bulk#Auto", 
            "version": "4.2.5"
        }
    }, 
    {
        "aggregation": "segment", 
        "analysis_id": "e2692740-f853-4dfa-99c2-3309f56cdfdb", 
        "basecall_1d": {
            "exit_status_dist": {
                "pass": 1
            }, 
            "qscore_dist_temp": [
                {
                    "count": 1, 
                    "mean_qscore": 13.0
                }
            ], 
            "qscore_sum_temp": {
                "count": 1, 
                "mean": 13.156586647033691, 
                "sum": 13.156586647033691
            }, 
            "read_len_events_sum_temp": 42150, 
            "seq_len_bases_dist_temp": [
                {
                    "count": 1, 
                    "length": 0.0
                }
            ], 
            "seq_len_bases_sum_temp": 1, 
            "seq_len_events_dist_temp": [
                {
                    "count": 1, 
                    "length": 42000.0
                }
            ], 
            "speed_bases_per_second_dist_temp": [
                {
                    "count": 1, 
                    "speed": 1.0
                }
            ], 
            "strand_median_pa": {
                "count": 1, 
                "mean": 70.73360443115234, 
                "sum": 70.73360443115234
            }, 
            "strand_sd_pa": {
                "count": 1, 
                "mean": 9.018534660339355, 
                "sum": 9.018534660339355
            }
        }, 
        "channel_count": 1, 
        "context_tags": {
            "barcoding_enabled": "0", 
            "experiment_duration_set": "1680", 
            "experiment_type": "genomic_dna", 
            "local_basecalling": "0", 
            "package": "bream4", 
            "package_version": "6.1.10", 
            "sample_frequency": "4000", 
            "sequencing_kit": "sqk-lsk110"
        }, 
        "latest_run_time": 71182.3359375, 
        "levels_sums": {
            "count": 1, 
            "mean": 186.80303955078125, 
            "open_pore_level_sum": 186.80303955078125
        }, 
        "opts": {
            "adapter_pt_range_scale": "5.200000", 
            "additional_lamp_context_bases": "2", 
            "align_ref": "", 
            "align_type": "auto", 
            "allow_inferior_barcodes": "0", 
            "as_cpu_threads_per_scaler": "2", 
            "as_gpu_runners_per_device": "2", 
            "as_model_file": "", 
            "as_num_scalers": "4", 
            "as_reads_per_runner": "32", 
            "bam_methylation_threshold": "5.000000", 
            "bam_out": "0", 
            "barcode_kits": "", 
            "barcode_nested_output_folder": "0", 
            "beam_cut": "100.000000", 
            "beam_width": "32", 
            "bed_file": "", 
            "builtin_scripts": "1", 
            "calib_detect": "0", 
            "calib_max_sequence_length": "3800", 
            "calib_min_coverage": "0.600000", 
            "calib_min_sequence_length": "3000", 
            "calib_reference": "lambda_3.6kb.fasta", 
            "chunk_size": "2000", 
            "chunks_per_caller": "10000", 
            "chunks_per_runner": "208", 
            "client_id": "-1", 
            "compress_fastq": "1", 
            "cpu_threads_per_caller": "4", 
            "detect_adapter": "0", 
            "detect_barcodes": "0", 
            "detect_mid_strand_adapter": "0", 
            "detect_mid_strand_barcodes": "0", 
            "detect_primer": "0", 
            "device": "cuda:0", 
            "disable_pings": "0", 
            "disable_qscore_filtering": "0", 
            "dmean_threshold": "1.000000", 
            "dmean_win_size": "2", 
            "do_read_splitting": "0", 
            "duplex_window_size_max": "1000", 
            "duplex_window_size_min": "200", 
            "end_gap1": "40", 
            "end_gap2": "40", 
            "extend_gap1": "40", 
            "extend_gap2": "160", 
            "fast5_out": "0", 
            "flowcell": "", 
            "front_window_size": "150", 
            "gpu_runners_per_device": "12", 
            "high_priority_threshold": "10", 
            "index": "0", 
            "input_file_list": "", 
            "int8_mode": "0", 
            "jump_threshold": "1.000000", 
            "kernel_path": "", 
            "kit": "", 
            "lamp_kit": "", 
            "log_speed_frequency": "0", 
            "max_queued_reads": "2000", 
            "max_read_split_depth": "2", 
            "max_search_len": "1000", 
            "medium_priority_threshold": "4", 
            "min_length_lamp_context": "30", 
            "min_length_lamp_target": "70", 
            "min_qscore": "10.000000", 
            "min_score_adapter": "60.000000", 
            "min_score_adapter_mid": "50.000000", 
            "min_score_barcode_front": "60.000000", 
            "min_score_barcode_mask": "40.000000", 
            "min_score_barcode_mid": "50.000000", 
            "min_score_barcode_rear": "60.000000", 
            "min_score_lamp": "80.000000", 
            "min_score_lamp_mask": "50.000000", 
            "min_score_lamp_target": "50.000000", 
            "min_score_primer": "60.000000", 
            "min_score_read_splitting": "70.000000", 
            "model_file": "template_r9.4.1_450bps_sup.jsn", 
            "moves_out": "0", 
            "nested_output_folder": "0", 
            "noisiest_section_scaling_max_size": "8000", 
            "num_alignment_threads": "4", 
            "num_barcode_threads": "4", 
            "num_barcoding_buffers": "24", 
            "num_base_mod_threads": "2", 
            "num_callers": "4", 
            "num_extra_bases_trim": "0", 
            "num_mid_barcoding_buffers": "96", 
            "num_read_splitting_buffers": "16", 
            "num_read_splitting_threads": "4", 
            "num_reads_per_barcoding_buffer": "4", 
            "open_gap1": "40", 
            "open_gap2": "160", 
            "overlap": "100", 
            "override_scaling": "0", 
            "ping_segment_duration": "60", 
            "ping_url": "https://ping.oxfordnanoportal.com/basecall", 
            "post_out": "0", 
            "print_workflows": "0", 
            "progress_stats_frequency": "-1.000000", 
            "pt_median_offset": "2.500000", 
            "pt_minimum_read_start_index": "30", 
            "pt_required_adapter_drop": "30.000000", 
            "pt_scaling": "0", 
            "qscore_offset": "0.349800", 
            "qscore_scale": "0.972200", 
            "quiet": "0", 
            "read_batch_size": "4000", 
            "read_id_list": "", 
            "read_splitting_arrangement_files": "", 
            "read_splitting_score_matrix_filename": "", 
            "rear_window_size": "150", 
            "records_per_fastq": "0", 
            "recursive": "0", 
            "remora_models": "", 
            "require_barcodes_both_ends": "0", 
            "resume": "0", 
            "reverse_sequence": "0", 
            "sample_sheet": "", 
            "scaling_mad": "1.000000", 
            "scaling_med": "0.000000", 
            "score_matrix_filename": "5x5_mismatch_matrix.txt", 
            "start_gap1": "40", 
            "start_gap2": "40", 
            "stay_penalty": "1.000000", 
            "temp_bias": "1.000000", 
            "temp_weight": "1.000000", 
            "trace_categories_logs": "", 
            "trace_domains_config": "", 
            "trim_adapters": "0", 
            "trim_barcodes": "0", 
            "trim_min_events": "3", 
            "trim_primers": "0", 
            "trim_strategy": "dna", 
            "trim_threshold": "2.500000", 
            "u_substitution": "0", 
            "verbose_logs": "0"
        }, 
        "read_count": 1, 
        "reads_per_channel_dist": [
            {
                "channel": 82, 
                "count": 1
            }
        ], 
        "run_id": "6a94eb36563a9975b38e81c41534f6590c24867c", 
        "segment_duration": 60, 
        "segment_number": 20, 
        "segment_type": "guppy-acquisition", 
        "software": {
            "analysis": "1d_basecalling", 
            "name": "guppy-basecalling", 
            "version": "6.1.2+e0556ff"
        }, 
        "tracking_id": {
            "asic_id": "3506339664", 
            "asic_id_eeprom": "5924305", 
            "asic_temp": "35.086052", 
            "asic_version": "IA02D", 
            "auto_update": "0", 
            "auto_update_source": "https://mirror.oxfordnanoportal.com/software/MinKNOW/", 
            "bream_is_standard": "0", 
            "configuration_version": "4.2.11", 
            "device_id": "MN34875", 
            "device_type": "minion", 
            "distribution_status": "stable", 
            "distribution_version": "21.02.1", 
            "exp_script_name": "sequencing_MIN106_DNA:FLO-FLG001:SQK-LSK110", 
            "exp_script_purpose": "sequencing_run", 
            "exp_start_time": "2021-08-03T19:16:19Z", 
            "flongle_adapter_id": "FA-00993", 
            "flow_cell_id": "AHC239", 
            "flow_cell_product_code": "FLO-FLG001", 
            "guppy_version": "4.3.4+ecb2805", 
            "heatsink_temp": "35.003906", 
            "hostname": "mckennalab", 
            "installation_type": "nc", 
            "local_firmware_file": "1", 
            "msg_id": "e18b1546-a071-4afd-9f9a-0ee10798498e", 
            "operating_system": "ubuntu 20.04", 
            "protocol_group_id": "FEX_120", 
            "protocol_run_id": "bdcd3859-b6c1-400c-b4e6-2bc287733f75", 
            "protocol_start_time": "", 
            "protocols_version": "6.1.10", 
            "run_id": "6a94eb36563a9975b38e81c41534f6590c24867c", 
            "sample_id": "96_v2", 
            "time_stamp": "2022-06-11T02:04:21Z", 
            "usb_config": "MinION_fx3_1.1.1_ONT#MinION_fpga_1.1.0#bulk#Auto", 
            "version": "4.2.5"
        }
    }, 
    {
        "aggregation": "segment", 
        "analysis_id": "e2692740-f853-4dfa-99c2-3309f56cdfdb", 
        "basecall_1d": {
            "exit_status_dist": {
                "pass": 2
            }, 
            "qscore_dist_temp": [
                {
                    "count": 1, 
                    "mean_qscore": 10.0
                }, 
                {
                    "count": 1, 
                    "mean_qscore": 12.0
                }
            ], 
            "qscore_sum_temp": {
                "count": 2, 
                "mean": 11.276811599731445, 
                "sum": 22.55362319946289
            }, 
            "read_len_events_sum_temp": 53138, 
            "seq_len_bases_dist_temp": [
                {
                    "count": 2, 
                    "length": 0.0
                }
            ], 
            "seq_len_bases_sum_temp": 2, 
            "seq_len_events_dist_temp": [
                {
                    "count": 1, 
                    "length": 22000.0
                }, 
                {
                    "count": 1, 
                    "length": 30000.0
                }
            ], 
            "speed_bases_per_second_dist_temp": [
                {
                    "count": 2, 
                    "speed": 1.0
                }
            ], 
            "strand_median_pa": {
                "count": 2, 
                "mean": 64.367584228515625, 
                "sum": 128.73516845703125
            }, 
            "strand_sd_pa": {
                "count": 2, 
                "mean": 8.753284454345703, 
                "sum": 17.506568908691406
            }
        }, 
        "channel_count": 2, 
        "context_tags": {
            "barcoding_enabled": "0", 
            "experiment_duration_set": "1680", 
            "experiment_type": "genomic_dna", 
            "local_basecalling": "0", 
            "package": "bream4", 
            "package_version": "6.1.10", 
            "sample_frequency": "4000", 
            "sequencing_kit": "sqk-lsk110"
        }, 
        "latest_run_time": 74867.7109375, 
        "levels_sums": {
            "count": 2, 
            "mean": 193.033935546875, 
            "open_pore_level_sum": 386.06787109375
        }, 
        "opts": {
            "adapter_pt_range_scale": "5.200000", 
            "additional_lamp_context_bases": "2", 
            "align_ref": "", 
            "align_type": "auto", 
            "allow_inferior_barcodes": "0", 
            "as_cpu_threads_per_scaler": "2", 
            "as_gpu_runners_per_device": "2", 
            "as_model_file": "", 
            "as_num_scalers": "4", 
            "as_reads_per_runner": "32", 
            "bam_methylation_threshold": "5.000000", 
            "bam_out": "0", 
            "barcode_kits": "", 
            "barcode_nested_output_folder": "0", 
            "beam_cut": "100.000000", 
            "beam_width": "32", 
            "bed_file": "", 
            "builtin_scripts": "1", 
            "calib_detect": "0", 
            "calib_max_sequence_length": "3800", 
            "calib_min_coverage": "0.600000", 
            "calib_min_sequence_length": "3000", 
            "calib_reference": "lambda_3.6kb.fasta", 
            "chunk_size": "2000", 
            "chunks_per_caller": "10000", 
            "chunks_per_runner": "208", 
            "client_id": "-1", 
            "compress_fastq": "1", 
            "cpu_threads_per_caller": "4", 
            "detect_adapter": "0", 
            "detect_barcodes": "0", 
            "detect_mid_strand_adapter": "0", 
            "detect_mid_strand_barcodes": "0", 
            "detect_primer": "0", 
            "device": "cuda:0", 
            "disable_pings": "0", 
            "disable_qscore_filtering": "0", 
            "dmean_threshold": "1.000000", 
            "dmean_win_size": "2", 
            "do_read_splitting": "0", 
            "duplex_window_size_max": "1000", 
            "duplex_window_size_min": "200", 
            "end_gap1": "40", 
            "end_gap2": "40", 
            "extend_gap1": "40", 
            "extend_gap2": "160", 
            "fast5_out": "0", 
            "flowcell": "", 
            "front_window_size": "150", 
            "gpu_runners_per_device": "12", 
            "high_priority_threshold": "10", 
            "index": "0", 
            "input_file_list": "", 
            "int8_mode": "0", 
            "jump_threshold": "1.000000", 
            "kernel_path": "", 
            "kit": "", 
            "lamp_kit": "", 
            "log_speed_frequency": "0", 
            "max_queued_reads": "2000", 
            "max_read_split_depth": "2", 
            "max_search_len": "1000", 
            "medium_priority_threshold": "4", 
            "min_length_lamp_context": "30", 
            "min_length_lamp_target": "70", 
            "min_qscore": "10.000000", 
            "min_score_adapter": "60.000000", 
            "min_score_adapter_mid": "50.000000", 
            "min_score_barcode_front": "60.000000", 
            "min_score_barcode_mask": "40.000000", 
            "min_score_barcode_mid": "50.000000", 
            "min_score_barcode_rear": "60.000000", 
            "min_score_lamp": "80.000000", 
            "min_score_lamp_mask": "50.000000", 
            "min_score_lamp_target": "50.000000", 
            "min_score_primer": "60.000000", 
            "min_score_read_splitting": "70.000000", 
            "model_file": "template_r9.4.1_450bps_sup.jsn", 
            "moves_out": "0", 
            "nested_output_folder": "0", 
            "noisiest_section_scaling_max_size": "8000", 
            "num_alignment_threads": "4", 
            "num_barcode_threads": "4", 
            "num_barcoding_buffers": "24", 
            "num_base_mod_threads": "2", 
            "num_callers": "4", 
            "num_extra_bases_trim": "0", 
            "num_mid_barcoding_buffers": "96", 
            "num_read_splitting_buffers": "16", 
            "num_read_splitting_threads": "4", 
            "num_reads_per_barcoding_buffer": "4", 
            "open_gap1": "40", 
            "open_gap2": "160", 
            "overlap": "100", 
            "override_scaling": "0", 
            "ping_segment_duration": "60", 
            "ping_url": "https://ping.oxfordnanoportal.com/basecall", 
            "post_out": "0", 
            "print_workflows": "0", 
            "progress_stats_frequency": "-1.000000", 
            "pt_median_offset": "2.500000", 
            "pt_minimum_read_start_index": "30", 
            "pt_required_adapter_drop": "30.000000", 
            "pt_scaling": "0", 
            "qscore_offset": "0.349800", 
            "qscore_scale": "0.972200", 
            "quiet": "0", 
            "read_batch_size": "4000", 
            "read_id_list": "", 
            "read_splitting_arrangement_files": "", 
            "read_splitting_score_matrix_filename": "", 
            "rear_window_size": "150", 
            "records_per_fastq": "0", 
            "recursive": "0", 
            "remora_models": "", 
            "require_barcodes_both_ends": "0", 
            "resume": "0", 
            "reverse_sequence": "0", 
            "sample_sheet": "", 
            "scaling_mad": "1.000000", 
            "scaling_med": "0.000000", 
            "score_matrix_filename": "5x5_mismatch_matrix.txt", 
            "start_gap1": "40", 
            "start_gap2": "40", 
            "stay_penalty": "1.000000", 
            "temp_bias": "1.000000", 
            "temp_weight": "1.000000", 
            "trace_categories_logs": "", 
            "trace_domains_config": "", 
            "trim_adapters": "0", 
            "trim_barcodes": "0", 
            "trim_min_events": "3", 
            "trim_primers": "0", 
            "trim_strategy": "dna", 
            "trim_threshold": "2.500000", 
            "u_substitution": "0", 
            "verbose_logs": "0"
        }, 
        "read_count": 2, 
        "reads_per_channel_dist": [
            {
                "channel": 76, 
                "count": 1
            }, 
            {
                "channel": 102, 
                "count": 1
            }
        ], 
        "run_id": "6a94eb36563a9975b38e81c41534f6590c24867c", 
        "segment_duration": 60, 
        "segment_number": 21, 
        "segment_type": "guppy-acquisition", 
        "software": {
            "analysis": "1d_basecalling", 
            "name": "guppy-basecalling", 
            "version": "6.1.2+e0556ff"
        }, 
        "tracking_id": {
            "asic_id": "3506339664", 
            "asic_id_eeprom": "5924305", 
            "asic_temp": "35.086052", 
            "asic_version": "IA02D", 
            "auto_update": "0", 
            "auto_update_source": "https://mirror.oxfordnanoportal.com/software/MinKNOW/", 
            "bream_is_standard": "0", 
            "configuration_version": "4.2.11", 
            "device_id": "MN34875", 
            "device_type": "minion", 
            "distribution_status": "stable", 
            "distribution_version": "21.02.1", 
            "exp_script_name": "sequencing_MIN106_DNA:FLO-FLG001:SQK-LSK110", 
            "exp_script_purpose": "sequencing_run", 
            "exp_start_time": "2021-08-03T19:16:19Z", 
            "flongle_adapter_id": "FA-00993", 
            "flow_cell_id": "AHC239", 
            "flow_cell_product_code": "FLO-FLG001", 
            "guppy_version": "4.3.4+ecb2805", 
            "heatsink_temp": "35.003906", 
            "hostname": "mckennalab", 
            "installation_type": "nc", 
            "local_firmware_file": "1", 
            "msg_id": "a3945ef7-0802-429d-949c-168e8a8e1885", 
            "operating_system": "ubuntu 20.04", 
            "protocol_group_id": "FEX_120", 
            "protocol_run_id": "bdcd3859-b6c1-400c-b4e6-2bc287733f75", 
            "protocol_start_time": "", 
            "protocols_version": "6.1.10", 
            "run_id": "6a94eb36563a9975b38e81c41534f6590c24867c", 
            "sample_id": "96_v2", 
            "time_stamp": "2022-06-11T02:04:21Z", 
            "usb_config": "MinION_fx3_1.1.1_ONT#MinION_fpga_1.1.0#bulk#Auto", 
            "version": "4.2.5"
        }
    }, 
    {
        "aggregation": "segment", 
        "analysis_id": "e2692740-f853-4dfa-99c2-3309f56cdfdb", 
        "basecall_1d": {
            "exit_status_dist": {
                "pass": 3
            }, 
            "qscore_dist_temp": [
                {
                    "count": 1, 
                    "mean_qscore": 11.5
                }, 
                {
                    "count": 2, 
                    "mean_qscore": 13.0
                }
            ], 
            "qscore_sum_temp": {
                "count": 3, 
                "mean": 12.72543716430664, 
                "sum": 38.17631149291992
            }, 
            "read_len_events_sum_temp": 98071, 
            "seq_len_bases_dist_temp": [
                {
                    "count": 3, 
                    "length": 0.0
                }
            ], 
            "seq_len_bases_sum_temp": 3, 
            "seq_len_events_dist_temp": [
                {
                    "count": 1, 
                    "length": 31000.0
                }, 
                {
                    "count": 2, 
                    "length": 33000.0
                }
            ], 
            "speed_bases_per_second_dist_temp": [
                {
                    "count": 3, 
                    "speed": 1.0
                }
            ], 
            "strand_median_pa": {
                "count": 3, 
                "mean": 71.26410675048828, 
                "sum": 213.79232788085938
            }, 
            "strand_sd_pa": {
                "count": 3, 
                "mean": 9.490092277526855, 
                "sum": 28.47027587890625
            }
        }, 
        "channel_count": 3, 
        "context_tags": {
            "barcoding_enabled": "0", 
            "experiment_duration_set": "1680", 
            "experiment_type": "genomic_dna", 
            "local_basecalling": "0", 
            "package": "bream4", 
            "package_version": "6.1.10", 
            "sample_frequency": "4000", 
            "sequencing_kit": "sqk-lsk110"
        }, 
        "latest_run_time": 85584.2265625, 
        "levels_sums": {
            "count": 3, 
            "mean": 206.9370574951172, 
            "open_pore_level_sum": 620.8111572265625
        }, 
        "opts": {
            "adapter_pt_range_scale": "5.200000", 
            "additional_lamp_context_bases": "2", 
            "align_ref": "", 
            "align_type": "auto", 
            "allow_inferior_barcodes": "0", 
            "as_cpu_threads_per_scaler": "2", 
            "as_gpu_runners_per_device": "2", 
            "as_model_file": "", 
            "as_num_scalers": "4", 
            "as_reads_per_runner": "32", 
            "bam_methylation_threshold": "5.000000", 
            "bam_out": "0", 
            "barcode_kits": "", 
            "barcode_nested_output_folder": "0", 
            "beam_cut": "100.000000", 
            "beam_width": "32", 
            "bed_file": "", 
            "builtin_scripts": "1", 
            "calib_detect": "0", 
            "calib_max_sequence_length": "3800", 
            "calib_min_coverage": "0.600000", 
            "calib_min_sequence_length": "3000", 
            "calib_reference": "lambda_3.6kb.fasta", 
            "chunk_size": "2000", 
            "chunks_per_caller": "10000", 
            "chunks_per_runner": "208", 
            "client_id": "-1", 
            "compress_fastq": "1", 
            "cpu_threads_per_caller": "4", 
            "detect_adapter": "0", 
            "detect_barcodes": "0", 
            "detect_mid_strand_adapter": "0", 
            "detect_mid_strand_barcodes": "0", 
            "detect_primer": "0", 
            "device": "cuda:0", 
            "disable_pings": "0", 
            "disable_qscore_filtering": "0", 
            "dmean_threshold": "1.000000", 
            "dmean_win_size": "2", 
            "do_read_splitting": "0", 
            "duplex_window_size_max": "1000", 
            "duplex_window_size_min": "200", 
            "end_gap1": "40", 
            "end_gap2": "40", 
            "extend_gap1": "40", 
            "extend_gap2": "160", 
            "fast5_out": "0", 
            "flowcell": "", 
            "front_window_size": "150", 
            "gpu_runners_per_device": "12", 
            "high_priority_threshold": "10", 
            "index": "0", 
            "input_file_list": "", 
            "int8_mode": "0", 
            "jump_threshold": "1.000000", 
            "kernel_path": "", 
            "kit": "", 
            "lamp_kit": "", 
            "log_speed_frequency": "0", 
            "max_queued_reads": "2000", 
            "max_read_split_depth": "2", 
            "max_search_len": "1000", 
            "medium_priority_threshold": "4", 
            "min_length_lamp_context": "30", 
            "min_length_lamp_target": "70", 
            "min_qscore": "10.000000", 
            "min_score_adapter": "60.000000", 
            "min_score_adapter_mid": "50.000000", 
            "min_score_barcode_front": "60.000000", 
            "min_score_barcode_mask": "40.000000", 
            "min_score_barcode_mid": "50.000000", 
            "min_score_barcode_rear": "60.000000", 
            "min_score_lamp": "80.000000", 
            "min_score_lamp_mask": "50.000000", 
            "min_score_lamp_target": "50.000000", 
            "min_score_primer": "60.000000", 
            "min_score_read_splitting": "70.000000", 
            "model_file": "template_r9.4.1_450bps_sup.jsn", 
            "moves_out": "0", 
            "nested_output_folder": "0", 
            "noisiest_section_scaling_max_size": "8000", 
            "num_alignment_threads": "4", 
            "num_barcode_threads": "4", 
            "num_barcoding_buffers": "24", 
            "num_base_mod_threads": "2", 
            "num_callers": "4", 
            "num_extra_bases_trim": "0", 
            "num_mid_barcoding_buffers": "96", 
            "num_read_splitting_buffers": "16", 
            "num_read_splitting_threads": "4", 
            "num_reads_per_barcoding_buffer": "4", 
            "open_gap1": "40", 
            "open_gap2": "160", 
            "overlap": "100", 
            "override_scaling": "0", 
            "ping_segment_duration": "60", 
            "ping_url": "https://ping.oxfordnanoportal.com/basecall", 
            "post_out": "0", 
            "print_workflows": "0", 
            "progress_stats_frequency": "-1.000000", 
            "pt_median_offset": "2.500000", 
            "pt_minimum_read_start_index": "30", 
            "pt_required_adapter_drop": "30.000000", 
            "pt_scaling": "0", 
            "qscore_offset": "0.349800", 
            "qscore_scale": "0.972200", 
            "quiet": "0", 
            "read_batch_size": "4000", 
            "read_id_list": "", 
            "read_splitting_arrangement_files": "", 
            "read_splitting_score_matrix_filename": "", 
            "rear_window_size": "150", 
            "records_per_fastq": "0", 
            "recursive": "0", 
            "remora_models": "", 
            "require_barcodes_both_ends": "0", 
            "resume": "0", 
            "reverse_sequence": "0", 
            "sample_sheet": "", 
            "scaling_mad": "1.000000", 
            "scaling_med": "0.000000", 
            "score_matrix_filename": "5x5_mismatch_matrix.txt", 
            "start_gap1": "40", 
            "start_gap2": "40", 
            "stay_penalty": "1.000000", 
            "temp_bias": "1.000000", 
            "temp_weight": "1.000000", 
            "trace_categories_logs": "", 
            "trace_domains_config": "", 
            "trim_adapters": "0", 
            "trim_barcodes": "0", 
            "trim_min_events": "3", 
            "trim_primers": "0", 
            "trim_strategy": "dna", 
            "trim_threshold": "2.500000", 
            "u_substitution": "0", 
            "verbose_logs": "0"
        }, 
        "read_count": 3, 
        "reads_per_channel_dist": [
            {
                "channel": 28, 
                "count": 1
            }, 
            {
                "channel": 82, 
                "count": 1
            }, 
            {
                "channel": 104, 
                "count": 1
            }
        ], 
        "run_id": "6a94eb36563a9975b38e81c41534f6590c24867c", 
        "segment_duration": 60, 
        "segment_number": 24, 
        "segment_type": "guppy-acquisition", 
        "software": {
            "analysis": "1d_basecalling", 
            "name": "guppy-basecalling", 
            "version": "6.1.2+e0556ff"
        }, 
        "tracking_id": {
            "asic_id": "3506339664", 
            "asic_id_eeprom": "5924305", 
            "asic_temp": "35.086052", 
            "asic_version": "IA02D", 
            "auto_update": "0", 
            "auto_update_source": "https://mirror.oxfordnanoportal.com/software/MinKNOW/", 
            "bream_is_standard": "0", 
            "configuration_version": "4.2.11", 
            "device_id": "MN34875", 
            "device_type": "minion", 
            "distribution_status": "stable", 
            "distribution_version": "21.02.1", 
            "exp_script_name": "sequencing_MIN106_DNA:FLO-FLG001:SQK-LSK110", 
            "exp_script_purpose": "sequencing_run", 
            "exp_start_time": "2021-08-03T19:16:19Z", 
            "flongle_adapter_id": "FA-00993", 
            "flow_cell_id": "AHC239", 
            "flow_cell_product_code": "FLO-FLG001", 
            "guppy_version": "4.3.4+ecb2805", 
            "heatsink_temp": "35.003906", 
            "hostname": "mckennalab", 
            "installation_type": "nc", 
            "local_firmware_file": "1", 
            "msg_id": "3fae6c83-c542-45cc-b66e-5a4c8dd93bf3", 
            "operating_system": "ubuntu 20.04", 
            "protocol_group_id": "FEX_120", 
            "protocol_run_id": "bdcd3859-b6c1-400c-b4e6-2bc287733f75", 
            "protocol_start_time": "", 
            "protocols_version": "6.1.10", 
            "run_id": "6a94eb36563a9975b38e81c41534f6590c24867c", 
            "sample_id": "96_v2", 
            "time_stamp": "2022-06-11T02:04:21Z", 
            "usb_config": "MinION_fx3_1.1.1_ONT#MinION_fpga_1.1.0#bulk#Auto", 
            "version": "4.2.5"
        }
    }, 
    {
        "aggregation": "segment", 
        "analysis_id": "e2692740-f853-4dfa-99c2-3309f56cdfdb", 
        "basecall_1d": {
            "exit_status_dist": {
                "pass": 1
            }, 
            "qscore_dist_temp": [
                {
                    "count": 1, 
                    "mean_qscore": 12.5
                }
            ], 
            "qscore_sum_temp": {
                "count": 1, 
                "mean": 12.603055000305176, 
                "sum": 12.603055000305176
            }, 
            "read_len_events_sum_temp": 30714, 
            "seq_len_bases_dist_temp": [
                {
                    "count": 1, 
                    "length": 0.0
                }
            ], 
            "seq_len_bases_sum_temp": 1, 
            "seq_len_events_dist_temp": [
                {
                    "count": 1, 
                    "length": 30000.0
                }
            ], 
            "speed_bases_per_second_dist_temp": [
                {
                    "count": 1, 
                    "speed": 1.0
                }
            ], 
            "strand_median_pa": {
                "count": 1, 
                "mean": 68.43476104736328, 
                "sum": 68.43476104736328
            }, 
            "strand_sd_pa": {
                "count": 1, 
                "mean": 9.018534660339355, 
                "sum": 9.018534660339355
            }
        }, 
        "channel_count": 1, 
        "context_tags": {
            "barcoding_enabled": "0", 
            "experiment_duration_set": "1680", 
            "experiment_type": "genomic_dna", 
            "local_basecalling": "0", 
            "package": "bream4", 
            "package_version": "6.1.10", 
            "sample_frequency": "4000", 
            "sequencing_kit": "sqk-lsk110"
        }, 
        "latest_run_time": 90691.9453125, 
        "levels_sums": {
            "count": 1, 
            "mean": 197.220458984375, 
            "open_pore_level_sum": 197.220458984375
        }, 
        "opts": {
            "adapter_pt_range_scale": "5.200000", 
            "additional_lamp_context_bases": "2", 
            "align_ref": "", 
            "align_type": "auto", 
            "allow_inferior_barcodes": "0", 
            "as_cpu_threads_per_scaler": "2", 
            "as_gpu_runners_per_device": "2", 
            "as_model_file": "", 
            "as_num_scalers": "4", 
            "as_reads_per_runner": "32", 
            "bam_methylation_threshold": "5.000000", 
            "bam_out": "0", 
            "barcode_kits": "", 
            "barcode_nested_output_folder": "0", 
            "beam_cut": "100.000000", 
            "beam_width": "32", 
            "bed_file": "", 
            "builtin_scripts": "1", 
            "calib_detect": "0", 
            "calib_max_sequence_length": "3800", 
            "calib_min_coverage": "0.600000", 
            "calib_min_sequence_length": "3000", 
            "calib_reference": "lambda_3.6kb.fasta", 
            "chunk_size": "2000", 
            "chunks_per_caller": "10000", 
            "chunks_per_runner": "208", 
            "client_id": "-1", 
            "compress_fastq": "1", 
            "cpu_threads_per_caller": "4", 
            "detect_adapter": "0", 
            "detect_barcodes": "0", 
            "detect_mid_strand_adapter": "0", 
            "detect_mid_strand_barcodes": "0", 
            "detect_primer": "0", 
            "device": "cuda:0", 
            "disable_pings": "0", 
            "disable_qscore_filtering": "0", 
            "dmean_threshold": "1.000000", 
            "dmean_win_size": "2", 
            "do_read_splitting": "0", 
            "duplex_window_size_max": "1000", 
            "duplex_window_size_min": "200", 
            "end_gap1": "40", 
            "end_gap2": "40", 
            "extend_gap1": "40", 
            "extend_gap2": "160", 
            "fast5_out": "0", 
            "flowcell": "", 
            "front_window_size": "150", 
            "gpu_runners_per_device": "12", 
            "high_priority_threshold": "10", 
            "index": "0", 
            "input_file_list": "", 
            "int8_mode": "0", 
            "jump_threshold": "1.000000", 
            "kernel_path": "", 
            "kit": "", 
            "lamp_kit": "", 
            "log_speed_frequency": "0", 
            "max_queued_reads": "2000", 
            "max_read_split_depth": "2", 
            "max_search_len": "1000", 
            "medium_priority_threshold": "4", 
            "min_length_lamp_context": "30", 
            "min_length_lamp_target": "70", 
            "min_qscore": "10.000000", 
            "min_score_adapter": "60.000000", 
            "min_score_adapter_mid": "50.000000", 
            "min_score_barcode_front": "60.000000", 
            "min_score_barcode_mask": "40.000000", 
            "min_score_barcode_mid": "50.000000", 
            "min_score_barcode_rear": "60.000000", 
            "min_score_lamp": "80.000000", 
            "min_score_lamp_mask": "50.000000", 
            "min_score_lamp_target": "50.000000", 
            "min_score_primer": "60.000000", 
            "min_score_read_splitting": "70.000000", 
            "model_file": "template_r9.4.1_450bps_sup.jsn", 
            "moves_out": "0", 
            "nested_output_folder": "0", 
            "noisiest_section_scaling_max_size": "8000", 
            "num_alignment_threads": "4", 
            "num_barcode_threads": "4", 
            "num_barcoding_buffers": "24", 
            "num_base_mod_threads": "2", 
            "num_callers": "4", 
            "num_extra_bases_trim": "0", 
            "num_mid_barcoding_buffers": "96", 
            "num_read_splitting_buffers": "16", 
            "num_read_splitting_threads": "4", 
            "num_reads_per_barcoding_buffer": "4", 
            "open_gap1": "40", 
            "open_gap2": "160", 
            "overlap": "100", 
            "override_scaling": "0", 
            "ping_segment_duration": "60", 
            "ping_url": "https://ping.oxfordnanoportal.com/basecall", 
            "post_out": "0", 
            "print_workflows": "0", 
            "progress_stats_frequency": "-1.000000", 
            "pt_median_offset": "2.500000", 
            "pt_minimum_read_start_index": "30", 
            "pt_required_adapter_drop": "30.000000", 
            "pt_scaling": "0", 
            "qscore_offset": "0.349800", 
            "qscore_scale": "0.972200", 
            "quiet": "0", 
            "read_batch_size": "4000", 
            "read_id_list": "", 
            "read_splitting_arrangement_files": "", 
            "read_splitting_score_matrix_filename": "", 
            "rear_window_size": "150", 
            "records_per_fastq": "0", 
            "recursive": "0", 
            "remora_models": "", 
            "require_barcodes_both_ends": "0", 
            "resume": "0", 
            "reverse_sequence": "0", 
            "sample_sheet": "", 
            "scaling_mad": "1.000000", 
            "scaling_med": "0.000000", 
            "score_matrix_filename": "5x5_mismatch_matrix.txt", 
            "start_gap1": "40", 
            "start_gap2": "40", 
            "stay_penalty": "1.000000", 
            "temp_bias": "1.000000", 
            "temp_weight": "1.000000", 
            "trace_categories_logs": "", 
            "trace_domains_config": "", 
            "trim_adapters": "0", 
            "trim_barcodes": "0", 
            "trim_min_events": "3", 
            "trim_primers": "0", 
            "trim_strategy": "dna", 
            "trim_threshold": "2.500000", 
            "u_substitution": "0", 
            "verbose_logs": "0"
        }, 
        "read_count": 1, 
        "reads_per_channel_dist": [
            {
                "channel": 43, 
                "count": 1
            }
        ], 
        "run_id": "6a94eb36563a9975b38e81c41534f6590c24867c", 
        "segment_duration": 60, 
        "segment_number": 26, 
        "segment_type": "guppy-acquisition", 
        "software": {
            "analysis": "1d_basecalling", 
            "name": "guppy-basecalling", 
            "version": "6.1.2+e0556ff"
        }, 
        "tracking_id": {
            "asic_id": "3506339664", 
            "asic_id_eeprom": "5924305", 
            "asic_temp": "35.086052", 
            "asic_version": "IA02D", 
            "auto_update": "0", 
            "auto_update_source": "https://mirror.oxfordnanoportal.com/software/MinKNOW/", 
            "bream_is_standard": "0", 
            "configuration_version": "4.2.11", 
            "device_id": "MN34875", 
            "device_type": "minion", 
            "distribution_status": "stable", 
            "distribution_version": "21.02.1", 
            "exp_script_name": "sequencing_MIN106_DNA:FLO-FLG001:SQK-LSK110", 
            "exp_script_purpose": "sequencing_run", 
            "exp_start_time": "2021-08-03T19:16:19Z", 
            "flongle_adapter_id": "FA-00993", 
            "flow_cell_id": "AHC239", 
            "flow_cell_product_code": "FLO-FLG001", 
            "guppy_version": "4.3.4+ecb2805", 
            "heatsink_temp": "35.003906", 
            "hostname": "mckennalab", 
            "installation_type": "nc", 
            "local_firmware_file": "1", 
            "msg_id": "ec5080db-cba6-423c-8dc5-ba25267933ba", 
            "operating_system": "ubuntu 20.04", 
            "protocol_group_id": "FEX_120", 
            "protocol_run_id": "bdcd3859-b6c1-400c-b4e6-2bc287733f75", 
            "protocol_start_time": "", 
            "protocols_version": "6.1.10", 
            "run_id": "6a94eb36563a9975b38e81c41534f6590c24867c", 
            "sample_id": "96_v2", 
            "time_stamp": "2022-06-11T02:04:21Z", 
            "usb_config": "MinION_fx3_1.1.1_ONT#MinION_fpga_1.1.0#bulk#Auto", 
            "version": "4.2.5"
        }
    }, 
    {
        "aggregation": "cumulative", 
        "analysis_id": "e2692740-f853-4dfa-99c2-3309f56cdfdb", 
        "basecall_1d": {
            "exit_status_dist": {
                "pass": 70
            }, 
            "qscore_dist_temp": [
                {
                    "count": 3, 
                    "mean_qscore": 10.0
                }, 
                {
                    "count": 4, 
                    "mean_qscore": 10.5
                }, 
                {
                    "count": 4, 
                    "mean_qscore": 11.0
                }, 
                {
                    "count": 8, 
                    "mean_qscore": 11.5
                }, 
                {
                    "count": 5, 
                    "mean_qscore": 12.0
                }, 
                {
                    "count": 13, 
                    "mean_qscore": 12.5
                }, 
                {
                    "count": 9, 
                    "mean_qscore": 13.0
                }, 
                {
                    "count": 11, 
                    "mean_qscore": 13.5
                }, 
                {
                    "count": 7, 
                    "mean_qscore": 14.0
                }, 
                {
                    "count": 5, 
                    "mean_qscore": 14.5
                }, 
                {
                    "count": 1, 
                    "mean_qscore": 15.0
                }
            ], 
            "qscore_sum_temp": {
                "count": 70, 
                "mean": 12.859704971313477, 
                "sum": 900.1793212890625
            }, 
            "read_len_events_sum_temp": 1831067, 
            "seq_len_bases_dist_temp": [
                {
                    "count": 70, 
                    "length": 0.0
                }
            ], 
            "seq_len_bases_sum_temp": 70, 
            "seq_len_events_dist_temp": [
                {
                    "count": 3, 
                    "length": 16000.0
                }, 
                {
                    "count": 3, 
                    "length": 17000.0
                }, 
                {
                    "count": 4, 
                    "length": 18000.0
                }, 
                {
                    "count": 1, 
                    "length": 19000.0
                }, 
                {
                    "count": 3, 
                    "length": 20000.0
                }, 
                {
                    "count": 2, 
                    "length": 21000.0
                }, 
                {
                    "count": 3, 
                    "length": 22000.0
                }, 
                {
                    "count": 9, 
                    "length": 23000.0
                }, 
                {
                    "count": 5, 
                    "length": 24000.0
                }, 
                {
                    "count": 7, 
                    "length": 25000.0
                }, 
                {
                    "count": 5, 
                    "length": 26000.0
                }, 
                {
                    "count": 5, 
                    "length": 27000.0
                }, 
                {
                    "count": 2, 
                    "length": 28000.0
                }, 
                {
                    "count": 1, 
                    "length": 29000.0
                }, 
                {
                    "count": 4, 
                    "length": 30000.0
                }, 
                {
                    "count": 3, 
                    "length": 31000.0
                }, 
                {
                    "count": 3, 
                    "length": 33000.0
                }, 
                {
                    "count": 1, 
                    "length": 34000.0
                }, 
                {
                    "count": 1, 
                    "length": 35000.0
                }, 
                {
                    "count": 2, 
                    "length": 38000.0
                }, 
                {
                    "count": 1, 
                    "length": 40000.0
                }, 
                {
                    "count": 1, 
                    "length": 42000.0
                }, 
                {
                    "count": 1, 
                    "length": 47000.0
                }
            ], 
            "speed_bases_per_second_dist_temp": [
                {
                    "count": 70, 
                    "speed": 1.0
                }
            ], 
            "strand_median_pa": {
                "count": 70, 
                "mean": 72.01439666748047, 
                "sum": 5041.0078125
            }, 
            "strand_sd_pa": {
                "count": 70, 
                "mean": 9.329257011413574, 
                "sum": 653.0479736328125
            }
        }, 
        "channel_count": 42, 
        "context_tags": {
            "barcoding_enabled": "0", 
            "experiment_duration_set": "1680", 
            "experiment_type": "genomic_dna", 
            "local_basecalling": "0", 
            "package": "bream4", 
            "package_version": "6.1.10", 
            "sample_frequency": "4000", 
            "sequencing_kit": "sqk-lsk110"
        }, 
        "latest_run_time": 90691.9453125, 
        "levels_sums": {
            "count": 70, 
            "mean": 193.89675903320312, 
            "open_pore_level_sum": 13572.7734375
        }, 
        "opts": {
            "adapter_pt_range_scale": "5.200000", 
            "additional_lamp_context_bases": "2", 
            "align_ref": "", 
            "align_type": "auto", 
            "allow_inferior_barcodes": "0", 
            "as_cpu_threads_per_scaler": "2", 
            "as_gpu_runners_per_device": "2", 
            "as_model_file": "", 
            "as_num_scalers": "4", 
            "as_reads_per_runner": "32", 
            "bam_methylation_threshold": "5.000000", 
            "bam_out": "0", 
            "barcode_kits": "", 
            "barcode_nested_output_folder": "0", 
            "beam_cut": "100.000000", 
            "beam_width": "32", 
            "bed_file": "", 
            "builtin_scripts": "1", 
            "calib_detect": "0", 
            "calib_max_sequence_length": "3800", 
            "calib_min_coverage": "0.600000", 
            "calib_min_sequence_length": "3000", 
            "calib_reference": "lambda_3.6kb.fasta", 
            "chunk_size": "2000", 
            "chunks_per_caller": "10000", 
            "chunks_per_runner": "208", 
            "client_id": "-1", 
            "compress_fastq": "1", 
            "cpu_threads_per_caller": "4", 
            "detect_adapter": "0", 
            "detect_barcodes": "0", 
            "detect_mid_strand_adapter": "0", 
            "detect_mid_strand_barcodes": "0", 
            "detect_primer": "0", 
            "device": "cuda:0", 
            "disable_pings": "0", 
            "disable_qscore_filtering": "0", 
            "dmean_threshold": "1.000000", 
            "dmean_win_size": "2", 
            "do_read_splitting": "0", 
            "duplex_window_size_max": "1000", 
            "duplex_window_size_min": "200", 
            "end_gap1": "40", 
            "end_gap2": "40", 
            "extend_gap1": "40", 
            "extend_gap2": "160", 
            "fast5_out": "0", 
            "flowcell": "", 
            "front_window_size": "150", 
            "gpu_runners_per_device": "12", 
            "high_priority_threshold": "10", 
            "index": "0", 
            "input_file_list": "", 
            "int8_mode": "0", 
            "jump_threshold": "1.000000", 
            "kernel_path": "", 
            "kit": "", 
            "lamp_kit": "", 
            "log_speed_frequency": "0", 
            "max_queued_reads": "2000", 
            "max_read_split_depth": "2", 
            "max_search_len": "1000", 
            "medium_priority_threshold": "4", 
            "min_length_lamp_context": "30", 
            "min_length_lamp_target": "70", 
            "min_qscore": "10.000000", 
            "min_score_adapter": "60.000000", 
            "min_score_adapter_mid": "50.000000", 
            "min_score_barcode_front": "60.000000", 
            "min_score_barcode_mask": "40.000000", 
            "min_score_barcode_mid": "50.000000", 
            "min_score_barcode_rear": "60.000000", 
            "min_score_lamp": "80.000000", 
            "min_score_lamp_mask": "50.000000", 
            "min_score_lamp_target": "50.000000", 
            "min_score_primer": "60.000000", 
            "min_score_read_splitting": "70.000000", 
            "model_file": "template_r9.4.1_450bps_sup.jsn", 
            "moves_out": "0", 
            "nested_output_folder": "0", 
            "noisiest_section_scaling_max_size": "8000", 
            "num_alignment_threads": "4", 
            "num_barcode_threads": "4", 
            "num_barcoding_buffers": "24", 
            "num_base_mod_threads": "2", 
            "num_callers": "4", 
            "num_extra_bases_trim": "0", 
            "num_mid_barcoding_buffers": "96", 
            "num_read_splitting_buffers": "16", 
            "num_read_splitting_threads": "4", 
            "num_reads_per_barcoding_buffer": "4", 
            "open_gap1": "40", 
            "open_gap2": "160", 
            "overlap": "100", 
            "override_scaling": "0", 
            "ping_segment_duration": "60", 
            "ping_url": "https://ping.oxfordnanoportal.com/basecall", 
            "post_out": "0", 
            "print_workflows": "0", 
            "progress_stats_frequency": "-1.000000", 
            "pt_median_offset": "2.500000", 
            "pt_minimum_read_start_index": "30", 
            "pt_required_adapter_drop": "30.000000", 
            "pt_scaling": "0", 
            "qscore_offset": "0.349800", 
            "qscore_scale": "0.972200", 
            "quiet": "0", 
            "read_batch_size": "4000", 
            "read_id_list": "", 
            "read_splitting_arrangement_files": "", 
            "read_splitting_score_matrix_filename": "", 
            "rear_window_size": "150", 
            "records_per_fastq": "0", 
            "recursive": "0", 
            "remora_models": "", 
            "require_barcodes_both_ends": "0", 
            "resume": "0", 
            "reverse_sequence": "0", 
            "sample_sheet": "", 
            "scaling_mad": "1.000000", 
            "scaling_med": "0.000000", 
            "score_matrix_filename": "5x5_mismatch_matrix.txt", 
            "start_gap1": "40", 
            "start_gap2": "40", 
            "stay_penalty": "1.000000", 
            "temp_bias": "1.000000", 
            "temp_weight": "1.000000", 
            "trace_categories_logs": "", 
            "trace_domains_config": "", 
            "trim_adapters": "0", 
            "trim_barcodes": "0", 
            "trim_min_events": "3", 
            "trim_primers": "0", 
            "trim_strategy": "dna", 
            "trim_threshold": "2.500000", 
            "u_substitution": "0", 
            "verbose_logs": "0"
        }, 
        "read_count": 70, 
        "reads_per_channel_dist": [
            {
                "channel": 9, 
                "count": 2
            }, 
            {
                "channel": 16, 
                "count": 3
            }, 
            {
                "channel": 17, 
                "count": 1
            }, 
            {
                "channel": 18, 
                "count": 1
            }, 
            {
                "channel": 23, 
                "count": 1
            }, 
            {
                "channel": 26, 
                "count": 2
            }, 
            {
                "channel": 27, 
                "count": 2
            }, 
            {
                "channel": 28, 
                "count": 3
            }, 
            {
                "channel": 29, 
                "count": 2
            }, 
            {
                "channel": 31, 
                "count": 3
            }, 
            {
                "channel": 32, 
                "count": 1
            }, 
            {
                "channel": 33, 
                "count": 1
            }, 
            {
                "channel": 34, 
                "count": 1
            }, 
            {
                "channel": 35, 
                "count": 1
            }, 
            {
                "channel": 36, 
                "count": 1
            }, 
            {
                "channel": 40, 
                "count": 1
            }, 
            {
                "channel": 42, 
                "count": 2
            }, 
            {
                "channel": 43, 
                "count": 1
            }, 
            {
                "channel": 49, 
                "count": 3
            }, 
            {
                "channel": 54, 
                "count": 3
            }, 
            {
                "channel": 55, 
                "count": 1
            }, 
            {
                "channel": 56, 
                "count": 1
            }, 
            {
                "channel": 62, 
                "count": 2
            }, 
            {
                "channel": 64, 
                "count": 2
            }, 
            {
                "channel": 66, 
                "count": 1
            }, 
            {
                "channel": 67, 
                "count": 1
            }, 
            {
                "channel": 71, 
                "count": 2
            }, 
            {
                "channel": 75, 
                "count": 1
            }, 
            {
                "channel": 76, 
                "count": 2
            }, 
            {
                "channel": 82, 
                "count": 3
            }, 
            {
                "channel": 83, 
                "count": 1
            }, 
            {
                "channel": 86, 
                "count": 1
            }, 
            {
                "channel": 89, 
                "count": 2
            }, 
            {
                "channel": 95, 
                "count": 2
            }, 
            {
                "channel": 99, 
                "count": 3
            }, 
            {
                "channel": 100, 
                "count": 1
            }, 
            {
                "channel": 102, 
                "count": 1
            }, 
            {
                "channel": 104, 
                "count": 1
            }, 
            {
                "channel": 109, 
                "count": 1
            }, 
            {
                "channel": 111, 
                "count": 4
            }, 
            {
                "channel": 112, 
                "count": 1
            }, 
            {
                "channel": 115, 
                "count": 1
            }
        ], 
        "run_id": "6a94eb36563a9975b38e81c41534f6590c24867c", 
        "segment_duration": 1560, 
        "segment_number": 1, 
        "segment_type": "guppy-acquisition", 
        "software": {
            "analysis": "1d_basecalling", 
            "name": "guppy-basecalling", 
            "version": "6.1.2+e0556ff"
        }, 
        "tracking_id": {
            "asic_id": "3506339664", 
            "asic_id_eeprom": "5924305", 
            "asic_temp": "35.086052", 
            "asic_version": "IA02D", 
            "auto_update": "0", 
            "auto_update_source": "https://mirror.oxfordnanoportal.com/software/MinKNOW/", 
            "bream_is_standard": "0", 
            "configuration_version": "4.2.11", 
            "device_id": "MN34875", 
            "device_type": "minion", 
            "distribution_status": "stable", 
            "distribution_version": "21.02.1", 
            "exp_script_name": "sequencing_MIN106_DNA:FLO-FLG001:SQK-LSK110", 
            "exp_script_purpose": "sequencing_run", 
            "exp_start_time": "2021-08-03T19:16:19Z", 
            "flongle_adapter_id": "FA-00993", 
            "flow_cell_id": "AHC239", 
            "flow_cell_product_code": "FLO-FLG001", 
            "guppy_version": "4.3.4+ecb2805", 
            "heatsink_temp": "35.003906", 
            "hostname": "mckennalab", 
            "installation_type": "nc", 
            "local_firmware_file": "1", 
            "msg_id": "6c473204-5ed3-492c-a1c5-068a8c372544", 
            "operating_system": "ubuntu 20.04", 
            "protocol_group_id": "FEX_120", 
            "protocol_run_id": "bdcd3859-b6c1-400c-b4e6-2bc287733f75", 
            "protocol_start_time": "", 
            "protocols_version": "6.1.10", 
            "run_id": "6a94eb36563a9975b38e81c41534f6590c24867c", 
            "sample_id": "96_v2", 
            "time_stamp": "2022-06-11T02:04:21Z", 
            "usb_config": "MinION_fx3_1.1.1_ONT#MinION_fpga_1.1.0#bulk#Auto", 
            "version": "4.2.5"
        }
    }
]