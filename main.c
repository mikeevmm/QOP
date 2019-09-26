#include "main.h"

int main() {
    printf("Hello world!");
}

/* 
impl<P> Gate<P> {
    fn new_from_matrix(matrix: Array2<Complex<f64>>, paramfn: Option<ReparamFn<P>>) -> Self {
        Gate {
            matrix,
            paramfn,
            id: GateId::Custom,
        }
    }

    fn new(gate: GateId) -> Self {
        let matrix = match (gate) {
            GateId::I => array![
                [Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)],
                [Complex::new(0.0, 0.0), Complex::new(1.0, 0.0)]
            ],
            GateId::X => array![
                [Complex::new(0.0, 0.0), Complex::new(1.0, 0.0)],
                [Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)]
            ],
            GateId::Y => array![
                [Complex::new(0.0, 0.0), Complex::new(0.0, -1.0)],
                [Complex::new(0.0, 1.0), Complex::new(1.0, 0.0)]
            ],
            GateId::Z => array![
                [Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)],
                [Complex::new(0.0, 0.0), Complex::new(-1.0, 0.0)]    
            ],
            GateId::H => 7.071067812*array![
                [Complex::new(1.0, 0.0), Complex::new(1.0, 0.0)],
                [Complex::new(1.0, 0.0), Complex::new(-1.0, 0.0)]
            ],
            GateId::SqrtX => 0.5*array![
                [Complex::new(1.0, 1.0), Complex::new(1.0, -1.0)],
                [Complex::new(1.0, -1.0), Complex::new(-1.0, 1.0)]
            ],
            GateId::T => 
            _ => unimplemented!()
        }
    }
}
 */