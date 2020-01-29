fn main() {
    println!("Hello, world!");
}

struct COF {
    a: Vec<f64>,
    b: Vec<f64>,
    n: usize,
}

// Assumes a is of size n*n and b is of size n
fn gauss(m: usize, cof: &mut COF) {
    println!("The value of m is: {}", m);
    let n = cof.b.len();

    for k in 1..n {
        let kp = k + 1;
        for i in kp..n+1 {
            let r = cof.a[i*n + k]/cof.a[k*n + k];
            for j in kp..n+1 {
                cof.a[i*n + j] -= r*cof.a[k*n + j];
            }
            for j in 1..m+1 {
                cof.b[i*n + j] -= r*cof.b[k*n + j];
            }
        }
    }

    for k in 1..m+1 {
        cof.b[n*n + k] /= cof.a[n*n + n];
        for i in (n-1..0).rev() {
            let ip = i + 1;
            for j in ip..n+1 {
                cof.b[i*n + k] -= cof.a[i*n + j] * cof.b[j*n + k];
            }
            cof.b[i*n + k] /= cof.a[i*n + i];
        }
    }

}