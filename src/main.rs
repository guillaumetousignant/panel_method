fn main() {
    println!("Hello, world!");
}

struct COF {
    A: Vec<f64>,
    B: Vec<f64>,
    N: usize,
}

// Assumes A is of size n*n and B is of size n
fn gauss(m: usize, cof: &mut COF) {
    println!("The value of m is: {}", m);
    let n = cof.B.len();

    for k in 1..n {
        let kp = k + 1;
        for i in kp..n+1 {
            let r = cof.A[i*n + k]/cof.A[k*n + k];
            for j in kp..n+1 {
                cof.A[i*n + j] -= r*cof.A[k*n + j];
            }
            for j in 1..m+1 {
                cof.B[i*n + j] -= r*cof.B[k*n + j];
            }
        }
    }

    for k in 1..m+1 {
        cof.B[n*n + k] /= cof.A[n*n + n];
        for i in (n-1..0).rev() {
            let ip = i + 1;
            for j in ip..n+1 {
                cof.B[i*n + k] -= cof.A[i*n + j] * cof.B[j*n + k];
            }
            cof.B[i*n + k] /= cof.A[i*n + i];
        }
    }

}