use crate::building_blocks::*;

#[cfg(feature = "parallel")]
use rayon::iter::*;

type Rez = (f64, Vec<f64>, Vec<f64>);

fn measure_mono_sound_window<T, U> (input: &[T]) -> Rez
where T: num::PrimInt + std::convert::TryFrom<usize> + Send + Sync, 
            U: morton_encoding::ValidKey<T> + Send + Sync,
{
    let coor_bits = (std::mem::size_of::<T>()*8) as u8;
    let key_bits = (std::mem::size_of::<U>()*8) as u8;
    let size = input.len();
    
    #[cfg(not(feature = "parallel"))]
    let clz_fun = get_clzs;
    #[cfg(feature = "parallel")]
    let clz_fun = get_clzs_par;

    let bloat_fn = |x| morton_encoding::bloat_custom::<T, U>(x, morton_encoding::nz(2));
    
    let tot = |x| -> T {T::from(x).unwrap()};
    let norm = |x| tot(x * (1<<coor_bits) / size);
    
    let encode_fn = |(a, &b)| (bloat_fn(norm(a))<<1) | bloat_fn(b);
    
    #[cfg(not(feature = "parallel"))]
    let buf = input.iter().enumerate();
    #[cfg(feature = "parallel")]
    let buf = input.par_iter().enumerate();
    
    let clzs = clz_fun(buf, encode_fn);
    let (tmp, lacun) = get_results_from_clzs(clzs, key_bits);
    finalise_results(tmp, lacun, size, key_bits, key_bits * 2)
}

fn measure_stereo_sound_window<T, U> (input: &[(T, T)]) -> Rez
where T: num::PrimInt + std::convert::TryFrom<usize> + Send + Sync, 
            U: morton_encoding::ValidKey<T> + Send + Sync,
{
    let coor_bits = (std::mem::size_of::<T>()*8) as usize;
    let key_bits = (std::mem::size_of::<U>()*8) as u8;
    let size = input.len();
    
    #[cfg(not(feature = "parallel"))]
    let clz_fun = get_clzs;
    #[cfg(feature = "parallel")]
    let clz_fun = get_clzs_par;

    let bloat_fn = |x| morton_encoding::bloat_custom::<T, U>(x, morton_encoding::nz(3));
    
    let tot = |x| -> T {T::from(x).unwrap()};
    let norm = |x| tot(x * (1<<coor_bits) / size);
    
    let encode_fn = |(a, &(b, c))| (bloat_fn(norm(a))<<2) | bloat_fn(b)<<1 | bloat_fn(c);
    
    #[cfg(not(feature = "parallel"))]
    let buf = input.iter().enumerate();
    #[cfg(feature = "parallel")]
    let buf = input.par_iter().enumerate();
    
    let clzs = clz_fun(buf, encode_fn);
    let (mut tmp, mut lacun) = get_results_from_clzs(clzs, key_bits);
    tmp.drain(0..coor_bits);
    lacun.drain(0..coor_bits);
    finalise_results(tmp, lacun, size, key_bits, key_bits * 3)
}


pub fn measure_sound_in_windows_mono (input: &[f64], window_size: usize) -> Vec<Rez> {
    
    if window_size < 256 {
        let frac_fun = |x: &[u8]| -> Rez {measure_mono_sound_window::<u8, u16>(x)};
        let norm_fun = |x: &f64| (*x * 256. / (window_size as f64)) as u8;
        let result: Vec<_> =    input.iter().map(norm_fun).collect();
        result.windows(window_size)
                .map(frac_fun)
                .collect()
    } else if window_size < (1<<16) {
        let frac_fun = |x: &[u16]| -> Rez {measure_mono_sound_window::<u16, u32>(x)};
        let norm_fun = |x: &f64| (*x * 65536. / (window_size as f64)) as u16;
        let result: Vec<_> =    input.iter().map(norm_fun).collect();
        result.windows(window_size)
                .map(frac_fun)
                .collect()
    } else {
        let frac_fun = |x: &[u32]| -> Rez {measure_mono_sound_window::<u32, u64>(x)};
        let norm_fun = |x: &f64| (*x * ((1u64<<32) as f64) / (window_size as f64)) as u32;
        let result: Vec<_> =    input.iter().map(norm_fun).collect();
        result.windows(window_size)
                .map(frac_fun)
                .collect()
    }
}


pub fn measure_sound_in_windows_stereo (input: &[f64], window_size: usize) -> Vec<Rez> {
    
    if window_size < 256 {
        let _frac_fun = |x: &[(u8, u8)]| -> Rez {measure_stereo_sound_window::<u8, u32>(x)};
    } drop((input, window_size)); unimplemented!()} /*
        let norm_fun = |x: &f64| (*x * 256. / (window_size as f64)) as u8;
        let result: Vec<_> =    input.iter().map(norm_fun).collect();
        result.chunks(2)
                .windows(window_size)
                .map(frac_fun)
                .collect()
    } else if window_size < (1<<16) {
        let frac_fun = |x: &[(u16, u16)]| -> Rez {measure_stereo_sound_window::<u16, u64>(x)};
        let norm_fun = |x: &f64| (*x * 65536. / (window_size as f64)) as u16;
        let result: Vec<_> =    input.iter().map(norm_fun).collect();
        result.chunks(2)
                .windows(window_size)
                .map(frac_fun)
                .collect()
    } else {
        let frac_fun = |x: &[(u32, u32)]| -> Rez {measure_stereo_sound_window::<u32, u128>(x)};
        let norm_fun = |x: &f64| (*x * ((1u64<<32) as f64) / (window_size as f64)) as u32;
        let result: Vec<_> =    input.iter().map(norm_fun).collect();
        result.chunks(2)
                .windows(window_size)
                .map(frac_fun)
                .collect()
    }
}*/
