pub struct BooleanFunction {
    truth_table: Vec<u8>,
    n: usize
}

impl BooleanFunction {

    fn new(truth_table: Vec<u8>) -> Self {
        let len = truth_table.len();
        assert!(len.is_power_of_two(), "Truth table length must be a power of two.");
        let n: usize = (usize::BITS - 1 - len.leading_zeros()) as usize;
        BooleanFunction { truth_table, n }
    }

    fn from_truth_table(truth_table: Vec<u8>) -> Self {
        Self::new(truth_table)
    }

    fn algebraic_immunity(&self) -> usize{
        crate::ai::AlgebraicImmunity::algebraic_immunity(self.truth_table.clone(), self.n)
    }
    
}