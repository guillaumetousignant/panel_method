fn main() {
    println!("Hello, world!");
}

struct BOD {
    ndtot: usize,
    x: Vec<f64>,
    y: Vec<f64>,
    x_mid: Vec<f64>,
    y_mid: Vec<f64>,
    costhe: Vec<f64>,
    sinthe: Vec<f64>,
}

struct COF {
    a: Vec<f64>,
    b: Vec<f64>,
    n: usize,
}

fn coef(sin_alpha: f64, cos_alpha: f64, bod: &mut BOD, cof: &mut COF) {

}   

// Assumes a is of size n*n and b is of size n
fn gauss(m: usize, cof: &mut COF) {
    let n = cof.b.len();

    for k in 0..n {
        let kp = k + 1;
        for i in kp..n+1 {
            let r = cof.a[i*n + k]/cof.a[k*n + k];
            for j in kp..n+1 {
                cof.a[i*n + j] -= r*cof.a[k*n + j];
            }
            for j in 0..m+1 {
                cof.b[i*n + j] -= r*cof.b[k*n + j];
            }
        }
    }

    for k in 0..m+1 {
        cof.b[n*n + k] /= cof.a[n*n + n];
        for i in (0..n).rev() {
            let ip = i + 1;
            for j in ip..n+1 {
                cof.b[i*n + k] -= cof.a[i*n + j] * cof.b[j*n + k];
            }
            cof.b[i*n + k] /= cof.a[i*n + i];
        }
    }

}