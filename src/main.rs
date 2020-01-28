fn main() {
    println!("Hello, world!");
}

// Assumes A is of size n*n and B is of size n
fn gauss(M: usize, A: &mut Vec<f64>, B: &mut Vec<f64>, n: isize) {
    println!("The value of M is: {}", M);
    let n = B.len();

    for k in 1..n {
        let kp = k + 1;
        for i in kp..n+1 {
            let r = A[i*n + k]/A[k*n + k];
            for j in kp..n+1 {
                A[i*n + j] -= r*A[k*n + j];
            }
            for j in 1..M+1 {
                B[i*n + j] -= r*B[k*n + j];
            }
        }
    }

    for k in 1..M+1 {
        B[n*n + k] /= A[n*n + n];
        for i in (n-1..0).rev() {
            let ip = i + 1;
            for j in ip..n+1 {
                B[i*n + k] -= A[i*n + j] * B[j*n + k];
            }
            B[i*n + k] /= A[i*n + i];
        }
    }

}