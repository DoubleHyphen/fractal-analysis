#[cfg(feature = "parallel")]
use rayon::prelude::*;


#[cfg(feature = "time_things")]
macro_rules! time {
    ($x: expr) => {{
        eprintln!("Measuring expression: {}", stringify!($x));
        let begin = std::time::Instant::now();
        let result = $x;
        eprintln!("Time elapsed: {:#?}\n", begin.elapsed());
        result
}};}

#[cfg(not(feature = "time_things"))]
macro_rules! time { ($x: expr) => {{$x}};}

#[cfg(feature = "parallel")]
pub fn get_clzs_par<H, Smp, Key, Set>(set: Set, get_key_from_sample: H) -> Vec<u8>
where
    H: Fn(Smp) -> Key + std::marker::Sync + std::marker::Send + Copy,
    Key: num_traits::int::PrimInt + std::marker::Sync + std::marker::Send,
    [Key] : rayon::prelude::ParallelSliceMut<Key>,
    Set : rayon::iter::IntoParallelIterator<Item = Smp> + std::clone::Clone,
    Smp: std::marker::Sync + std::marker::Send
{
    let mut tmp = time! (set
            .into_par_iter()
            .map(get_key_from_sample)
            .collect::<Vec<_>>());
        
        time!(    tmp.par_sort_unstable()   );
        
        time!(    tmp
        .par_windows(2)
        .map(|x| x[0] ^ x[1])
        .map(|x| x.leading_zeros() as u8)
        .collect()    )
}

pub fn get_clzs <H, Smp, Key, Set>(set: Set, get_key_from_sample: H) -> Vec<u8>
where
    H: Fn(Smp) -> Key,
    Key: num_traits::int::PrimInt,
    Set: IntoIterator<Item = Smp>,
{
        let mut tmp = time! (set
            .into_iter()
            .map(get_key_from_sample)
            .collect::<Vec<_>>());
        
        time!(    tmp.sort_unstable()   );
        
        time!(    tmp
        .windows(2)
        .map(|x| x[0] ^ x[1])
        .map(|x| x.leading_zeros() as u8)
        .collect()    )
}

pub fn get_inclination(input: &[f64]) -> f64
{
    let length = input.iter().count() as f64;
    let avy: f64 = input.iter().sum::<f64>()/length;
    let avx: f64 = length*(length-1.0)/(2.0 * length);
    let num_inc = |(i, x): (usize, f64)| -> f64 {(x-avy)*(i as f64 - avx)};
    let denom_inc = |(i, _): (usize, f64)| -> f64 {(i as f64 - avx)*(i as f64 - avx)};
    let num: f64 = input.iter().enumerate().map(|(a, b): (usize, &f64)| num_inc((a, *b))).sum();
    let denom: f64 = input.iter().enumerate().map(|(a, b): (usize, &f64)| denom_inc((a, *b))).sum();
    num/denom
}

pub fn get_results_from_clzs(input: Vec<u8>, key_bit_amt: u8) -> (Vec<u32>, Vec<u64>) {
    let mut s: Vec<u32> = vec![0; key_bit_amt as usize];
    let mut prevs: Vec<usize> = vec![0; key_bit_amt as usize];
    let mut squares: Vec<u64> = vec![0; key_bit_amt as usize];
    for (i, x) in input.into_iter().chain(vec![0].into_iter()).enumerate() {
        for b_i in (x as usize)..(key_bit_amt as usize) {
            s[b_i] += 1;
            squares[b_i] += (i - prevs[b_i]) as u64 * (i - prevs[b_i]) as u64;
            prevs[b_i] = i;
        }
    }
    (s, squares)
}

pub fn zbox_merge<H, Smp, Key, Set>(set: Set, get_key_from_sample: H, key_bits: u8, sample_size: usize, coor_bits: u8) -> (f64, Vec<f64>, Vec<f64>)
where
    H: Fn(Smp) -> Key,
    Key: num_traits::int::PrimInt,
    Set: IntoIterator<Item = Smp>,
{
    let clzs = get_clzs(set, get_key_from_sample);
    let (s, squares) = get_results_from_clzs(clzs, key_bits);
    finalise_results(s, squares, sample_size, coor_bits, key_bits)
}

#[cfg(feature = "parallel")]
pub fn zbox_merge_par<H, Smp, Key, Set>(set: Set, get_key_from_sample: H, key_bits: u8, sample_size: usize, coor_bits: u8) -> (f64, Vec<f64>, Vec<f64>)
where
    H: Fn(Smp) -> Key + std::marker::Sync + std::marker::Send + Copy,
    Key: num_traits::int::PrimInt + std::marker::Sync + std::marker::Send,
    [Key] : rayon::prelude::ParallelSliceMut<Key>,
    Set : rayon::iter::IntoParallelIterator<Item = Smp> + std::clone::Clone,
    Smp: std::marker::Sync + std::marker::Send
{
    let clzs = get_clzs_par(set, get_key_from_sample);
    let (s, squares) = get_results_from_clzs(clzs, key_bits);
    finalise_results(s, squares, sample_size, coor_bits, key_bits)
}
/*
pub fn get_results_from_clzs_functional(input: Vec<u8>, key_bit_amt: u8) -> (Vec<u32>, Vec<u64>) {
    let get_one_size = |x| {
        let tmp = input.iter()
            .chain(vec![0].iter())
            .enumerate()
            .filter(|(a, b)| b>=x);
        let size = tmp.len();
    };
    unimplemented!()
    let mut s: Vec<u32> = vec![0; key_bit_amt as usize];
    let mut prevs: Vec<usize> = vec![0; key_bit_amt as usize];
    let mut squares: Vec<u64> = vec![0; key_bit_amt as usize];
    for (i, x) in input.into_iter().chain(vec![0].into_iter()).enumerate() {
        for b_i in (x as usize)..(key_bit_amt as usize) {
            s[b_i] += 1;
            squares[b_i] += (i - prevs[b_i]) as u64 * (i - prevs[b_i]) as u64;
            prevs[b_i] = i;
        }
    }
    (s, squares)
}*/

pub fn finalise_results(s: Vec<u32>, squares: Vec<u64>, sample_size: usize, coor_bit_amt: u8, key_bit_amt: u8) -> (f64, Vec<f64>, Vec<f64>) {
    let step = (key_bit_amt/coor_bit_amt) as usize;
    let result_2 = s.iter().skip(step-1).step_by(step).map(|&x| f64::from(x).log2()).collect::<Vec<_>>();
    let result_3 = squares.into_iter().zip(s.into_iter()).skip(step-1).step_by(step).map(|(a, b)| (a as f64)*(b as f64)/(sample_size as f64 * sample_size as f64) - 1.0).collect::<Vec<_>>();
    let cap = (sample_size as f64).log2();
    let result_1_lim = result_2.iter().position(|x| *x>(0.9)*cap).unwrap_or(coor_bit_amt as usize);
    let result_1 = get_inclination(&result_2[0..result_1_lim]);
    (result_1, result_2, result_3)
}

