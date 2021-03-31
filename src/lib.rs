//! Welcome to the `fractal-analysis` crate! Here you will find necessary tools
//! for the analysis of fractal sets, based on the “Z-box merging” algorithm.
//!
//! # How to use
//! I have tried to create ergonomic functions for the following data-sets:
//! * [Images, colour or gray-scale](mod@img_funs)
//! * [Sounds](mod@snd_funs) (to be added later)
//! * [Electro-encephalograms](mod@eeg_19c_funs) (to be added later)
//! * [MRI images](mod@mri_funs) (to be added later)
//!
//! For each of those, you will find ready-made functions that can give you
//! useful results; just eg say `example_code_here` and you will get the result.
//! 
//! # Caveat
//! Currently (March 2021) the documentation and testing are still unfinished.
//! They are both provided on a best-effort basis; the user's understanding is
//! appreciated.
//!
//! # How to extend to other uses
//! Oh, how I _wish_ there was a way for me to design an API that said “Oh, just
//! do `zbox_merge(input)` for whatever input and it'll work!”.
//!
//! Sadly, it is no-where nearly as simple as that. The
//! basic library will get you a long way
//! along, but there is serious work you'll need to do before-hand. The simplest
//! way would be to just
//! [send the crate maintainer an e-mail](mailto:velona@ahiru.eu)
//! and ask directly, but in case s/he is unreachable, here is what you'll
//! need to do:
//!
//! ### Choose your dependent and independent dimensions
//! First, decide what counts as an “separate dimension” for your use-case.
//! Is each channel of your EEG a separate dimension that depends on the
//! time-stamp, or are they a sampling
//! of a two-dimensional scalp? What about the colour channels of your image?
//! That will give you a first idea of what to do, and denote the limits of the
//! result you will end up with.
//!
//! Of particular note: You will want the number of dimensions to be independent
//! of the measurement method. For instance: For the colour channels, we chose
//! each to be a separate dimension in part because all colour images have those
//! exact same 3 channels. But for the EEGs we chose to consider them as a 2-D
//! grid, or else the result would depend on the number of electrodes and
//! measurements between different setups would have no immediate way to be
//! compared. As a result, about half the mental effort was spent in deciding
//! exactly how to subdivide the electrode placement into quadrants!
//!
//! ### Count the number of bits you have available
//! For each coordinate, count the amount of values it can take. Take the
//! smallest of those amounts, find its base-2 logarithm,
//! and throw away all the decimal digits. The result
//! is the number of bits to which each of your coordinates will need to be
//! normalised.
//!
//! ### Choose your data-types
//! This is only important if you _really_ care about efficiency. The default
//! choices ought to be optimal for most cases, but under special circumstances
//! they _might_ not be.
//!
//! **Example 1:** If you have 5 keys of 12 bits each, each `Coor` will be an
//! `u16`, and as a result the default `Key` will be an `u128`. But they can
//! actually fit within an `u64`; the compiler just doesn't know that. In that
//! case, you might prefer to explicitly cast all keys to `u64`s.
//!
//! **Example 2:** If you need 48 bits for your key, and you're running the
//! code on a 16-bit architecture CPU, you might prefer to implement
//! custom-made 48 bit data-types (3*16) instead of letting the code use `u64`s
//! by default. The [`lindel`](lindel) crate permits that.
//!  
//! ### Create a “key-extraction” function
//! The key-extraction is comprised of two parts: Normalisation and
//! linearisation.
//!
//! The normalisation is already provided for you, and only
//! requires you to provide, for each co-ordinate, the min/max values and the
//! amount of bits to which they will be normalised. However, please bear in
//! mind that you will need to extract the independent coordinates together
//! with the dependent ones!
//!
//! As for the linearisation, two separate methods are already provided, thanks
//! to the [`lindel`](lindel) crate. Those are Morton encoding (“z-index”) and Hilbert
//! encoding (“Hilbert index”). If you care about speed and/or have
//! at least one independent coordinate, the z-index will serve you better.
//! If you have no independent coordinates, and you can afford the small
//! performance hit, a Hilbert index will easily allow you to examine almost
//! every possible subdivision of your data-set independently. Still, please
//! note the following:
//! 1. You will need to do that manually, as the program can't perform this
//! operation automatically.
//! 2. You could theoretically use a Hilbert index everywhere, but in
//! subdividing it you might end up omitting parts of the independent
//! coordinates whose values are too large.
//! 3. While a z-index is much quicker to compute than a Hilbert index, for
//! large data-sets the most expensive operation is the sorting, so perhaps the
//! difference won't be important for the `O(N)` part of the program.
//!
//! ### Choose a sorting algorithm
//! If you don't have a great concern for speed, you could choose the default
//! methods defined in `std` or `rayon`. If you need to squeeze every last bit
//! of performance, on the other hand, you might prefer to implement a radix
//! sort or something to that extent.
//!
//! ### And you're finished!
//! …finished with the preliminary work, that is. Now all that's left 
//! is to write the actual code: extract the keys, sort, extract the
//! amounts of common leading bits between each consecutive pair,
//! get the results from them, drain the useless bits, and use a
//! regression function to get the result. Here is some example code,
//! copied from the `img_funs` module:
//! ```should_not_compile
//!let width = input.width() as usize;
//!let height = input.height() as usize;
//!
//!let get_key_from_sample = |(x, y, rgb): (u32, u32, &image::Rgb<u8>)| -> u64 {
//!    let norm_x = create_normaliser::<u32, u8>(0, input.width()).unwrap();
//!    // The provided normaliser simplifies things slightly.
//!    let norm_y = create_normaliser::<u32, u8>(0, input.height()).unwrap();
//!    let arr = [norm_x(x), norm_y(y), rgb[0], rgb[1], rgb[2]];
//!    // Extract both independent and dependent coordinates
//!    morton_encoding::morton_encode(arr)
//!    // Linearise using `morton_encoding`
//!};
//!
//!let clzs = get_clzs(input.enumerate_pixels(), get_key_from_sample);
//!let (tmp, lacun) = get_results_from_clzs(clzs);
//!let tmp = drain_useless_bits::<_, 64, 40>(tmp);
//!let lacun = drain_useless_bits::<_, 64, 40>(lacun);
//!finalise_results(tmp, lacun, width*height, 8)
//! ```
//!
//! Unless of course you need more than 256 bits for each key, in which case
//! you'll need to change the rest of the functions so that they operate on
//! `u16`s instead.
//!
//! # Interpreting the results
//! ## Fractal dimension
//! _The first_ and most important _thing to understand_ with regards to the
//! fractal dimension is the limits inside which its
//! log-log plot will fall. Those create a quadrilateral, whose four edges
//! are defined as follows:
//! 1. Edge 1 is vertical, for `x` equal to the amount of bits you have
//! available. You can't subdivide your data-set more then the resolution that
//! each coordinate has.
//! 2. Edge 2 is horizontal, for `y` equal to the amount of samples you have
//! available. Subdivide all you might, you will never get more populated
//! boxes than the amount of samples (eg pixels) that exist in the data-set.
//! 3. Edge 3 is diagonal, and has a slope equal to the amount of independent
//! coordinates, if any.
//! 4. Edge 4 is diagonal as well; its slope is equal to the total amount
//! of coordinates available, independent as well as dependent. A data-set
//! constrained within `N` dimensions can never have a fractal dimension
//! above `N`.
//!
//! The most immediate result of these four constraints is that, especially if
//! the log-log plot rises abruptly at first, it will encounter a plateau as
//! soon as it reaches the total amount of samples, past which it will be
//! horizontal. As such, the horizontal part of the plot _must_ be excised
//! before calculating the fractal dimension, else it will be artificially low.
//!
//! _The second thing to understand_ is that, although an ideally fractal
//! shape's log-log plot will be linear, in practice data-sets will be
//! scale-dependent, leading to non-linearities in the log-log plot. In every
//! such instance we've found, the log-log plot is concave downward. Therefore,
//! the user has two choices:
//! 1. To take the simple linear regression of the log-log plot, and interpret
//! its slope as the fractal dimension. The scale-dependence, if any, may be
//! found from the mean square error between the line and the actual plot.
//! 2. To take a parabolic regression of the log-log plot, interpret the linear
//! parametre as the fractal dimension, and the second-order parametre as the
//! scale-dependence.
//!
//! Neither has been tried with any real rigour; the user is cordially invited
//! to share the results, if any.
//!
//! ## Lacunarity
//! The way the lacunarity was defined in the available literature, it appears
//! to be a measure that's different for each scale. It's measured here for the
//! user's convenience, but the interpretation of the results must necessarily
//! be left up to the user.
//!
//!  

#![cfg_attr(not(feature = "std"), no_std)]

use core::ops::{Add, BitXor, Div, Sub};
use num_traits::{CheckedShl, CheckedSub, PrimInt, Zero};
use arrayvec::ArrayVec;


pub mod img_funs;
//pub mod snd_funs;
//pub mod eeg_19c_funs;

pub use img_funs::*;
//pub use snd_funs::*;
//pub use eeg_19c_funs::*;

#[cfg(feature = "time_things")]
macro_rules! time {
    ($x: expr) => {{
        eprintln!("Measuring expression: {}", stringify!($x));
        let begin = std::time::Instant::now();
        let result = $x;
        eprintln!("Time elapsed: {:#?}\n", begin.elapsed());
        result
    }};
}

#[cfg(not(feature = "time_things"))]
macro_rules! time {
    ($x: expr) => {{
        $x
    }};
}

/// A convenience function that takes an arbitrary value, along with its minimum and maximum value, and normalises it to the limits of the coordinate data-type it's given.
/// ```rust
/// # use fractal_analysis::normalise;
/// let x = 1024u32;
/// let norm_x: u8 = normalise(x, 0, 65536);
/// assert_eq!(norm_x, 4);
/// ```
/// Please bear in mind:
/// 1. This function assumes that the minimum will be inclusive and the maximum exclusive.
/// 2. This function will silently return a zero if `maximum * (Coordinate::MAX + 1)` overflows, or if the value provided is outside `minimum..maximum`.
pub fn normalise<Input, Coordinate>(x: Input, minimum: Input, maximum: Input) -> Coordinate
where
    Input: Add<Output = Input>
        + Sub<Output = Input>
        + Div<Output = Input>
        + core::convert::TryInto<Coordinate>
        + CheckedShl
        + CheckedSub
        + Copy,
    Coordinate: Zero,
{
    let spread: Input = maximum - minimum;
    let coor_bits = core::mem::size_of::<Coordinate>() * 8;
    x.checked_sub(&minimum)
        .and_then(|x| x.checked_shl(coor_bits as u32))
        .map(|x| x / spread)
        .and_then(|x| x.try_into().ok())
        .unwrap_or(Coordinate::zero())
}

/// A convenience function that takes a minimum and maximum value, inclusive and exclusive respectively, and returns a function to normalise values to the limits of the coordinate data-type it's given.
/// ```rust
/// # use fractal_analysis::create_normaliser;
/// let norm_fn = create_normaliser::<u32, u8>(0u32, 8192).unwrap();
/// assert_eq!(norm_fn(31), 0);
/// assert_eq!(norm_fn(32), 1);
/// assert_eq!(norm_fn(250), 7);
/// ```
/// Please note:
/// 1. This function returns `None` if `maximum * (Coordinate::MAX + 1)` overflows.
/// 2. The output function will return a zero if given a value that's outside bounds.
pub fn create_normaliser<Input, Coordinate>(
    minimum: Input,
    maximum: Input,
) -> Option<impl Fn(Input) -> Coordinate>
where
    Input: Add<Output = Input>
        + Sub<Output = Input>
        + Div<Output = Input>
        + core::convert::TryInto<Coordinate>
        + CheckedShl
        + CheckedSub
        + Copy,
    Coordinate: Zero,
{
    let coor_bits = core::mem::size_of::<Coordinate>() * 8;
    let _ = maximum.checked_shl(coor_bits as u32)?;
    Some(move |x| normalise(x, minimum, maximum))
}

pub fn map_sampler<'a, H: 'a, Smp, Key: 'a, Set>(
    set: Set,
    get_key_from_sample: H,
) -> impl Iterator<Item = Key> + 'a
where
    H: Fn(Smp) -> Key,
    Key: PrimInt,
    Set: IntoIterator<Item = Smp> + 'a,
{
    set.into_iter().map(get_key_from_sample)
}

#[cfg(feature = "std")]
pub fn iterate_sorted_pairs<Key: Ord + Copy>(
    input: impl Iterator<Item = Key>,
) -> impl Iterator<Item = (Key, Key)> {
    let mut imp = input.collect::<Vec<_>>();
    imp.sort_unstable();
    (0..imp.len().saturating_sub(1)).map(move |i| (imp[i], imp[i + 1]))
}

#[cfg(not(feature = "std"))]
pub fn iterate_sorted_pairs<'a, Key: Ord + Copy>(
    input: impl Iterator<Item = Key>,
    buffer: &'a mut [Key],
) -> impl Iterator<Item = (Key, Key)> + 'a 
{
    buffer.iter_mut()
        .zip(input)
        .for_each(|(a, b)| *a = b);
    buffer.sort_unstable();
    (0..buffer.len().saturating_sub(1)).map(move |i| (buffer[i], buffer[i + 1]))
}

// TODO: Do this as a method to allow chaining

#[cfg(feature = "std")]
pub fn get_clzs<'a, H: 'a, Smp: 'a, Key: 'a, Set: 'a>(
    set: Set,
    get_key_from_sample: H,
) -> impl Iterator<Item = u8> + 'a
where
    H: Fn(Smp) -> Key,
    Key: PrimInt + BitXor<Output = Key>,
    Set: IntoIterator<Item = Smp>,
{
    let keys = time!(map_sampler(set, get_key_from_sample));
    let sorted_pairs_of_keys = time!(iterate_sorted_pairs(keys));
    time!(sorted_pairs_of_keys
        .map(|(a, b)| a ^ b)
        .map(|x| x.leading_zeros() as u8))
}

#[cfg(not(feature = "std"))]
pub fn get_clzs<'a, H: 'a, Smp: 'a, Key: 'a, Set: 'a>(
    set: Set,
    get_key_from_sample: H,
    buffer: &'a mut [Key]
) -> impl Iterator<Item = u8> + 'a
where
    H: Fn(Smp) -> Key,
    Key: PrimInt + BitXor<Output = Key>,
    Set: IntoIterator<Item = Smp>,
{
    let keys = time!(map_sampler(set, get_key_from_sample));
    let sorted_pairs_of_keys = time!(iterate_sorted_pairs(keys, buffer));
    time!(sorted_pairs_of_keys
        .map(|(a, b)| a ^ b)
        .map(|x| x.leading_zeros() as u8))
}

pub fn get_results_from_clzs<const KEY_BIT_AMT: usize>(
    input: impl Iterator<Item = u8>,
) -> (ArrayVec<u32, KEY_BIT_AMT>, ArrayVec<u64, KEY_BIT_AMT>) {
    let mut s = ArrayVec::from([0u32; KEY_BIT_AMT]);
    let mut prevs = ArrayVec::from([0usize; KEY_BIT_AMT]);
    let mut squares = ArrayVec::from([0u64; KEY_BIT_AMT]);
    for (i, x) in input.into_iter().chain(Some(0).into_iter()).enumerate() {
        for b_i in (x as usize)..(KEY_BIT_AMT) {
            s[b_i] += 1;
            squares[b_i] += (i - prevs[b_i]) as u64 * (i - prevs[b_i]) as u64;
            prevs[b_i] = i;
        }
    }
    (s, squares)
}

pub fn get_results_from_clzs_functional<const KEY_BIT_AMT: usize>(
    input: impl Iterator<Item = u8> + Clone,
) -> (ArrayVec<u32, KEY_BIT_AMT>, ArrayVec<usize, KEY_BIT_AMT>)
{
    let mut s = ArrayVec::from([0u32; KEY_BIT_AMT]);
    let mut squares = ArrayVec::from([0usize; KEY_BIT_AMT]);
    let smaller_clz_positions = |x| input.clone().enumerate().filter(move |(_, a)| *a <= x).map(|x| x.0);
    
    s.iter_mut().enumerate().for_each(|(i, x)| {
        *x = (smaller_clz_positions(i as u8)).count() as u32
    });
    
    squares.iter_mut().enumerate().for_each(|(i, a)| {
            let poss_1 = smaller_clz_positions(i as u8);
            let poss_2 = smaller_clz_positions(i as u8).chain(Some(0).into_iter()).skip(1);
            *a = poss_2.zip(poss_1).map(|x| x.0 - x.1).map(|x| x*x).sum();
        });
        
    (s, squares)
}

//const fn get_results_from_clzs = get_results_from_clzs_imperative;

pub fn get_inclination(input: &[f64]) -> f64 {
    let length = input.iter().count() as f64;
    let avy: f64 = input.iter().sum::<f64>() / length;
    let avx: f64 = length * (length - 1.0) / (2.0 * length);
    let num_inc = |(i, x): (usize, f64)| -> f64 { (x - avy) * (i as f64 - avx) };
    let denom_inc = |(i, _): (usize, f64)| -> f64 { (i as f64 - avx) * (i as f64 - avx) };
    let num: f64 = input
        .iter()
        .enumerate()
        .map(|(a, b): (usize, &f64)| num_inc((a, *b)))
        .sum();
    let denom: f64 = input
        .iter()
        .enumerate()
        .map(|(a, b): (usize, &f64)| denom_inc((a, *b)))
        .sum();
    num / denom
}

pub fn finalise_results<const KEY_BIT_AMT: usize>(
    s: ArrayVec<u32, KEY_BIT_AMT>,
    squares: ArrayVec<u64, KEY_BIT_AMT>,
    sample_size: usize,
    coor_bit_amt: u8,
) -> (f64, ArrayVec<f64, KEY_BIT_AMT>, ArrayVec<f64, KEY_BIT_AMT>) {
    let step = KEY_BIT_AMT / (coor_bit_amt as usize);
    let result_2 = s
        .iter()
        .skip(step - 1)
        .step_by(step)
        .map(|&x| f64::from(x).log2())
        .collect::<ArrayVec<_, KEY_BIT_AMT>>();
    let result_3 = squares
        .iter()
        .zip(s.iter())
        .skip(step - 1)
        .step_by(step)
        .map(|(&a, &b)| (a as f64) * (b as f64) / (sample_size as f64 * sample_size as f64) - 1.0)
        .collect::<ArrayVec<_, KEY_BIT_AMT>>();
    let cap = (sample_size as f64).log2();
    let result_1_lim = result_2
        .iter()
        .position(|x| *x > (0.9) * cap)
        .unwrap_or(coor_bit_amt as usize);
    let result_1 = get_inclination(&result_2[0..result_1_lim]);
    (result_1, result_2, result_3)
}

#[cfg(feature = "std")]
pub fn zbox_merge<H, Smp, Key, Set, const KEY_BIT_AMT: usize>(set: Set, get_key_from_sample: H, sample_size: usize, coor_bits: u8) -> (f64, ArrayVec<f64, KEY_BIT_AMT>, ArrayVec<f64, KEY_BIT_AMT>)
where
    H: Fn(Smp) -> Key,
    Key: PrimInt,
    Set: IntoIterator<Item = Smp>,
{
    let clzs = get_clzs(set, get_key_from_sample);
    let (s, squares) = get_results_from_clzs(clzs);
    finalise_results(s, squares, sample_size, coor_bits)
}

pub fn drain_useless_bits<T, const TOTAL_BITS: usize, const USEFUL_BITS: usize>(mut input: ArrayVec<T, TOTAL_BITS>) -> ArrayVec<T, USEFUL_BITS> {
    let useless_bits = TOTAL_BITS - USEFUL_BITS;
    input.drain(..useless_bits);
    input.into_iter().collect()
}






#[cfg(feature = "parallel")]
use rayon::prelude::*;

#[cfg(feature = "parallel")]
pub fn get_clzs_par<'a, H:'a, Smp:'a, Key:'a, Set:'a>(set: Set, get_key_from_sample: H) -> impl rayon::iter::ParallelIterator<Item = u8>
where
    H: Fn(Smp) -> Key + std::marker::Sync + std::marker::Send + Copy,
    Key: num_traits::int::PrimInt + std::marker::Sync + std::marker::Send,
    [Key] : rayon::prelude::ParallelSliceMut<Key>,
    Set : rayon::iter::IntoParallelIterator<Item = Smp>,
    Smp: std::marker::Sync + std::marker::Send
{
    let mut buffer = time! (set
            .into_par_iter()
            .map(get_key_from_sample)
            .collect::<Vec<_>>());

        time!(
            buffer.par_sort_unstable()
        );

        time!(
            (0..buffer.len().saturating_sub(1))
            .into_par_iter()
            .map(move |i| (buffer[i], buffer[i + 1]))
            .map(|x| x.0 ^ x.1)
            .map(|x| x.leading_zeros() as u8)
        )
}

#[cfg(feature = "parallel")]
pub fn get_results_from_clzs_parallel<const KEY_BIT_AMT: usize>(
    input: impl Iterator<Item = u8> + Clone + Sync,
) -> ([u32; KEY_BIT_AMT], [usize; KEY_BIT_AMT])
{
    let mut s = [0u32; KEY_BIT_AMT];
    let mut squares = [0usize; KEY_BIT_AMT];
    let smaller_clz_positions = |x| input.clone().enumerate().filter(move |(_, a)| *a <= x).map(|x| x.0);
    
    s.par_iter_mut().enumerate().for_each(|(i, x)| {
        *x = (smaller_clz_positions(i as u8)).count() as u32
    });
    
    squares.par_iter_mut().enumerate().for_each(|(i, a)| {
            let poss_1 = smaller_clz_positions(i as u8);
            let poss_2 = smaller_clz_positions(i as u8).chain(Some(0).into_iter()).skip(1);
            *a = poss_2.zip(poss_1).map(|x| x.0 - x.1).map(|x| x*x).sum();
        });
        
    (s, squares)
}

#[cfg(feature = "parallel")]
pub fn zbox_merge_par<H, Smp, Key, Set, const KEY_BIT_AMT: usize>(set: Set, get_key_from_sample: H, sample_size: usize, coor_bits: u8) -> (f64, ArrayVec<f64, KEY_BIT_AMT>, ArrayVec<f64, KEY_BIT_AMT>)
where
    H: Fn(Smp) -> Key + std::marker::Sync + std::marker::Send + Copy,
    Key: num_traits::int::PrimInt + std::marker::Sync + std::marker::Send,
    [Key] : rayon::prelude::ParallelSliceMut<Key>,
    Set : rayon::iter::IntoParallelIterator<Item = Smp>,
    Smp: std::marker::Sync + std::marker::Send
{
    let fnoo = get_results_from_clzs;
    let clzs = get_clzs_par(set, get_key_from_sample).collect::<Vec<_>>();
    let (s, squares) = fnoo(clzs.into_iter());
    finalise_results(s, squares, sample_size, coor_bits)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wtf() {
        assert!(true);
    }

    #[test]
    fn img_sorta() {
        use itertools::Itertools;
        use noise::utils::*;
        use noise::Fbm;
        use noise::MultiFractal;

        let mut dur = std::time::Duration::from_millis(0);
        for x in 0..10 {
            let clamp = |x, lower, upper| {
                if x < lower {
                    lower
                } else if x > upper {
                    upper
                } else {
                    x
                }
            };
            let size_x = 1024;
            let size_y = size_x;
            let fbm = Fbm::new().set_octaves(32).set_persistence(x as f64 / 10.);
            let pmb = PlaneMapBuilder::new(&fbm)
                .set_size(size_x as usize, size_y as usize)
                .build();
            let pix_span = (0..size_x as usize).cartesian_product(0..size_y as usize);
            let mut result = image::GrayImage::new(size_x as u32, size_y as u32);
            pix_span
                .clone()
                .map(|(xx, yy)| pmb.get_value(xx, yy))
                .map(|val| (val + 1.0) * 128.0)
                .map(|val| (clamp(val, 0.0, 255.99)) as u8)
                .zip(pix_span)
                .for_each(|(val, (xx, yy))| {
                    result.put_pixel(xx as u32, yy as u32, image::Luma([val]))
                });
            drop((pmb, fbm));

            let begin = std::time::Instant::now();
            let (hrm, _, _) = measure_lum_2d_u8_image(result);
            dur += begin.elapsed();
            dbg!(hrm);
        }
        println!("Time elapsed: {:?}", dur);

        use noise::Seedable;
        let mut dur = std::time::Duration::from_millis(0);
        for x in 0..10 {
            let clamp = |x, lower, upper| {
                if x < lower {
                    lower
                } else if x > upper {
                    upper
                } else {
                    x
                }
            };
            let size_x = 1 << 10;
            let size_y = size_x;
            let fbm_1 = Fbm::new()
                .set_octaves(32)
                .set_persistence(x as f64 / 10.)
                .set_seed(0);
            let fbm_2 = Fbm::new()
                .set_octaves(32)
                .set_persistence(x as f64 / 10.)
                .set_seed(1);
            let fbm_3 = Fbm::new()
                .set_octaves(32)
                .set_persistence(x as f64 / 10.)
                .set_seed(2);
            let pmb_1 = PlaneMapBuilder::new(&fbm_1)
                .set_size(size_x as usize, size_y as usize)
                .build();
            let pmb_2 = PlaneMapBuilder::new(&fbm_2)
                .set_size(size_x as usize, size_y as usize)
                .build();
            let pmb_3 = PlaneMapBuilder::new(&fbm_3)
                .set_size(size_x as usize, size_y as usize)
                .build();
            let pix_span = (0..size_x as usize).cartesian_product(0..size_y as usize);
            let mut result = image::RgbImage::new(size_x as u32, size_y as u32);
            let map_fn = |val: f64| (val + 1.0) * 128.0;
            let map_fn = |val: f64| clamp(map_fn(val), 0.0, 255.99) as u8;
            let map_fn = |x: [f64; 3]| [map_fn(x[0]), map_fn(x[1]), map_fn(x[2])];
            pix_span
                .clone()
                .map(|(xx, yy)| {
                    [
                        pmb_1.get_value(xx, yy),
                        pmb_2.get_value(xx, yy),
                        pmb_3.get_value(xx, yy),
                    ]
                })
                .map(map_fn)
                .zip(pix_span)
                .for_each(|(val, (xx, yy))| {
                    result.put_pixel(xx as u32, yy as u32, image::Rgb(val))
                });
            drop((pmb_1, fbm_1, pmb_2, fbm_2, pmb_3, fbm_3));

            let begin = std::time::Instant::now();
            let (hrm, _, _) = measure_rgb_2d_u8_image(result);
            dur += begin.elapsed();
            dbg!(hrm);
        }
        println!("Time elapsed: {:?}", dur);
    }
}
