use super::*;

#[cfg(feature = "parallel")]
use rayon::{iter::*, slice::*};

#[cfg(not(feature = "parallel"))]
pub fn measure_rgb_2d_u8_image (input: image::RgbImage) -> (f64, ArrayVec<f64, 40>, ArrayVec<f64, 40>) {
    let width = input.width() as usize;
    let height = input.height() as usize;

    let get_key_from_sample = |(x, y, rgb): (u32, u32, &image::Rgb<u8>)| -> u64 {
        let norm_x = create_normaliser::<u32, u8>(0, input.width()).unwrap();
        let norm_y = create_normaliser::<u32, u8>(0, input.height()).unwrap();
        let arr = [norm_x(x), norm_y(y), rgb[0], rgb[1], rgb[2]];
        morton_encoding::morton_encode(arr)
    };
    
    let clzs = get_clzs(input.enumerate_pixels(), get_key_from_sample);
    let (tmp, lacun) = get_results_from_clzs(clzs);
    let tmp = drain_useless_bits::<_, 64, 40>(tmp);
    let lacun = drain_useless_bits::<_, 64, 40>(lacun);
    finalise_results(tmp, lacun, width*height, 8)
}

#[cfg(not(feature = "parallel"))]
pub fn measure_lum_2d_u8_image (input: image::GrayImage) -> (f64, ArrayVec<f64, 24>, ArrayVec<f64, 24>) {
    let width = input.width() as usize;
    let height = input.height() as usize;

    let get_key_from_sample = |(x, y, lum): (u32, u32, &image::Luma<u8>)| -> u32 {
        let normalise_as_u8 = |q, max, min| ((q - min) * 256 / max) as u8;
        let norm_x = |x| normalise_as_u8(x, input.width(), 0);
        let norm_y = |y| normalise_as_u8(y, input.height(), 0);
        let arr = [norm_x(x), norm_y(y), lum[0]];
        morton_encoding::morton_encode(arr)
    };
    
    let clzs = get_clzs(input.enumerate_pixels(), get_key_from_sample);
    let (tmp, lacun) = get_results_from_clzs(clzs);
    let tmp = drain_useless_bits::<_, 32, 24>(tmp);
    let lacun = drain_useless_bits::<_, 32, 24>(lacun);
    finalise_results(tmp, lacun, width*height, 8)
}

#[cfg(feature = "parallel")]
pub fn measure_lum_2d_u8_image (input: image::GrayImage) -> (f64, ArrayVec<f64, 24>, ArrayVec<f64, 24>) {
    let width = input.width() as usize;
    let height = input.height() as usize;
    let buf = input.into_raw().into_par_iter();
    // Sadly, since `image` doesn't offer parallel iterators on its own,
    // we have to imitate this ourselves.
    let buf = buf.enumerate();
    
    let get_key_from_sample = |(comb, lum): (usize, u8)| -> u32 {
        let normalise_as_u8 = |q, max, min| ((q - min) * 256 / max) as u8;
        let norm_x = |x| normalise_as_u8(x, width, 0);
        let norm_y = |y| normalise_as_u8(y, height, 0);
        let x = comb % width;
        let y = comb / width;
        let x = norm_x(x);
        let y = norm_y(y);
        let arr = [x, y, lum];
        morton_encoding::morton_encode(arr)
    };
    
    let clzs = get_clzs_par(buf, get_key_from_sample).collect::<Vec<_>>();
    let (tmp, lacun) = get_results_from_clzs(clzs.into_iter());
    let tmp = drain_useless_bits::<_, 32, 24>(tmp);
    let lacun = drain_useless_bits::<_, 32, 24>(lacun);
    finalise_results(tmp, lacun, width*height, 8)
}

 
#[cfg(feature = "parallel")]
pub fn measure_rgb_2d_u8_image (input: image::RgbImage) -> (f64, ArrayVec<f64, 40>, ArrayVec<f64, 40>) {
    
    
// y as usize*self.width as usize + x as usize
    let width = input.width() as usize;
    let height = input.height() as usize;
    let buf = input.into_raw();
    // Sadly, since `image` doesn't offer parallel iterators on its own,
    // we have to imitate this ourselves.
    let buf = buf.par_chunks(3).enumerate();
    
    let get_key_from_sample = |(comb, rgb): (usize, &[u8])| -> u64 {
        let normalise_as_u8 = |q, max, min| ((q - min) * 256 / max) as u8;
        let norm_x = |x| normalise_as_u8(x, width, 0);
        let norm_y = |y| normalise_as_u8(y, height, 0);
        let x = norm_x(comb % width);
        let y = norm_y(comb / width);
        let arr = [x, y, rgb[0], rgb[1], rgb[2]];
        morton_encoding::morton_encode(arr)
    };
    
    let clzs = get_clzs_par(buf, get_key_from_sample).collect::<Vec<_>>();
    let (tmp, lacun) = get_results_from_clzs(clzs.into_iter());
    let tmp = drain_useless_bits::<_, 64, 40>(tmp);
    let lacun = drain_useless_bits::<_, 64, 40>(lacun);
    finalise_results(tmp, lacun, width*height, 8)
}
