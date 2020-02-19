use crate::building_blocks::*;

use std::fs::File;
use std::io::prelude::*;
use std::io;
use itertools::Itertools;

/** Tools for fractal analysis of ElectroEncephaloGrams (EEGs)

Still in a very rough shape, will update later. Nonetheless, it can be used as-is for EEGs of 19 channels.*/


const fn log_4_floor_times_2(chan_amt: usize) -> u8 {
    let lz = (chan_amt as u8).leading_zeros() as u8;
    7 - (lz | 1)
}

/**The maximal amount of {samples per fractal measurement} that make sense, given the number of electrodes.

In this case, the EEG samples can be sub-divided much more finely than the number of coordinates can, so the amount of EEG electrodes is what defines the maximal amount of samples we can take for each fractal measurement. */
const fn window_size(chan_amt: usize) -> usize {
    let immediately_lower_power_of_four = 1<<log_4_floor_times_2(chan_amt);
    immediately_lower_power_of_four
}

const fn normalisation_quotient(chan_amt: usize) -> usize {
    let tmp = (window_size(chan_amt).trailing_zeros())>>1;
    1<<tmp
}

const CHAN_AMT: usize = 19;
      
/** Reads `window_size` samples from a given EEG.

Assumes the following:
* Each EEG channel is in a different file
* Said files are named identically, apart from a serial number
* The numbers run from `1` to `chan_amt`.*/
pub fn read_window_not_normalised(folder: &str, filename: &str, extension: &str, window_size: usize, offset: usize, chan_amt: usize) -> Vec<Vec<f64>> {
    let cont_tmp = || -> io::Result<Vec<Vec<f64>>> {
        let mut entire_contents: Vec<Vec<f64>> = Vec::new();
        for serial_number in 1..=chan_amt {
            let filename = format!("{}{}{}{}", folder, filename, serial_number, extension);
            let mut f = File::open(filename)?;
        
            let mut sootooroo = String::new();
            f.read_to_string(&mut sootooroo)?;
            let partial_contents = sootooroo
                            .split_whitespace()
                            .filter_map(|x| x.parse::<f64>().ok())
                            .skip(offset)
                            .take(window_size)
                            .collect_vec();
            entire_contents.push(partial_contents);
        }
        Ok(entire_contents)
    };
    let mut result = cont_tmp().unwrap_or(Vec::new());
    let min_length = result.iter().map(|x| x.len()).min();
    if let Some(min_length) = min_length {
        for subvec in &mut result {
            subvec.truncate(min_length);
        }
    }
    result
}

// Heuristic assignment of x - y coordinates to the 19 electrodes of the EEGs we had to measure.
const COORS: [(u8, u8); CHAN_AMT] = [(0,1), (0,2), (1,0), (0,2), (2,1), (1,2), (3,1), (2,2), (3,1), (3,3), (0,0), (0,3), (2,0), (1,3), (3,0), (2,3), (1,1), (2,2), (3,2)];

/** Extracts a sub-window that's short enough (depending on the amount of channels of the EEG) and then normalises every sample depending on that specific sub-window length. For 19 electrodes, this leads to 16 samples per channel, normalised between 0 and 3 inclusive.
*/
pub fn extract_normalised_subwindow(input: &[Vec<f64>], offset: usize) -> Vec<Vec<u8>> {
    let channel_amt = input.len();
    let window_size = window_size(channel_amt);
    let sub_range = offset..(offset + window_size);
    let sub_window: Vec<_> = input[sub_range].to_owned();
    let norm_quot = (normalisation_quotient(channel_amt) as f64) - 0.0000001;
    
    let get_min_max = |x: &Vec<f64>| {
        let (&a, &b) = x.iter()
        .minmax()
        .into_option()
        .expect("Can't find min/max of empty iterator");
        (a, b)
    };
    
    sub_window
        .into_iter()
        .map(|x: Vec<f64>| {
            let (min, max) = get_min_max(&x);
            let norm_fun = |xx: f64| ((norm_quot*(xx-min)/(max-min)).floor()) as u8;
            x.into_iter().map(norm_fun).collect()
        }).collect()
}

type Rez = (f64, Vec<f64>, Vec<f64>);

pub fn get_one_result (x: Vec<Vec<u8>>) -> Rez {
    let channel_amt = x.len();
    let sample_amt = x[0].len();
    let sample_size = channel_amt * sample_amt;
    let coor_bits = log_4_floor_times_2(channel_amt)>>1;
    let norm_quot = (normalisation_quotient(channel_amt)) as f64;
    let mut normalised_samples_all_coors: Vec<[u8; 4]> = Vec::with_capacity(channel_amt * sample_amt);
    for (channel_id, channel_time_series) in x.into_iter().enumerate() {
        for (timestamp, normalised_sample) in channel_time_series.into_iter().enumerate() {
            let (x_coor, y_coor) = COORS[channel_id];
            let normalised_timestamp = (norm_quot * timestamp as f64 / sample_amt as f64) as u8;
            normalised_samples_all_coors.push([x_coor, y_coor, normalised_timestamp, normalised_sample]);
        }
    }
    let clzs = get_clzs(normalised_samples_all_coors, morton_encoding::morton_encode_u8_4d);
    let (s, squares) = get_results_from_clzs(clzs, coor_bits * 4);
    
    let result_2 = s.iter().map(|&x| f64::from(x).log2()).collect::<Vec<_>>();
    let result_3 = squares
            .into_iter()
            .zip(s.into_iter())
            .map(|(a, b)| ((a as f64)*(b as f64)/(sample_size as f64 * sample_size as f64)) - 1.0)
            .collect::<Vec<_>>();
    let result_1 = get_inclination(&result_2);
    (result_1, result_2, result_3)
}

pub fn get_result_vector_partial (folder: &str, filename: &str, extension: &str, window_size: usize, offset: usize, chan_amt: usize) -> Vec<Rez> {
    let window_nonn = read_window_not_normalised(folder, filename, extension, window_size, offset, chan_amt);
    let actual_window_size = window_nonn[0].len();
    let one_timestamp_result = |x| -> Rez {
        get_one_result(extract_normalised_subwindow(&window_nonn, x))
    };
    let limit = actual_window_size - window_size + 1;
    (0..limit).map(one_timestamp_result).collect()
}

pub fn get_result_vector_everything(folder: &str, filename: &str, extension: &str, chan_amt: usize) -> Vec<Rez> {
    get_result_vector_partial(folder, filename, extension, 0, 0usize.wrapping_sub(1), chan_amt)
}
