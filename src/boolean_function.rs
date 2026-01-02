use crate::ai::{AlgebraicImmunity, RestrictedAlgebraicImmunity};

pub struct BooleanFunction {
    t_table: Vec<u8>,
    n: usize,
}

impl BooleanFunction {
    /// Creates a Boolean function from its truth table.
    ///
    /// # Panics
    /// Panics if the truth table length is not a power of two.
    pub fn new(truth_table: Vec<u8>) -> Self {
        let len = truth_table.len();
        assert!(
            len.is_power_of_two(),
            "Truth table length must be a power of two."
        );
        let n: usize = (usize::BITS - 1 - len.leading_zeros()) as usize;
        BooleanFunction {
            t_table: truth_table,
            n: n,
        }
    }

    /// Alias for [`BooleanFunction::new`]
    ///
    /// # Arguments
    ///
    /// * `truth_table` - truth table of the Boolean function
    ///
    /// # Returns
    ///
    /// Boolean function object corresponding to `truth_table`.
    ///
    /// # Panics
    ///
    /// Panics if truth table lenght is not a power of 2.
    ///
    /// # Example
    ///
    /// ```
    /// use algebraic_immunity::boolean_function::BooleanFunction;
    /// let truth_table = vec![0,1,1,0];
    /// let boolean_function = BooleanFunction::from_truth_table(truth_table.clone());
    /// assert_eq!(boolean_function.truth_table(), truth_table);
    /// ```
    pub fn from_truth_table(truth_table: Vec<u8>) -> Self {
        Self::new(truth_table)
    }

    /// Get the Boolean function truth table from the Boolean function.
    ///
    /// # Returns
    ///
    /// Truth table of the Boolean function.
    ///
    /// # Example
    /// ```
    /// use algebraic_immunity::boolean_function::BooleanFunction;
    /// let t_table = vec![0,1,1,0];
    /// let bf = BooleanFunction::new(t_table.clone());
    /// assert_eq!(bf.truth_table(), t_table);
    /// ```
    pub fn truth_table(&self) -> Vec<u8> {
        self.t_table.clone()
    }

    /// Computes the algebraic immunity of the Boolean function.
    ///
    /// # Returns
    ///
    /// The algebraic immunity of the Boolean function, denoted `AI(f)`.
    ///
    /// # Examples
    ///
    /// ```
    /// use algebraic_immunity::boolean_function::BooleanFunction;
    ///
    /// let t_table = vec![0, 1, 1, 0];
    /// let bf = BooleanFunction::new(t_table);
    /// assert_eq!(bf.algebraic_immunity(), 1);
    /// ```
    pub fn algebraic_immunity(&self) -> usize {
        AlgebraicImmunity::algebraic_immunity(self.t_table.clone(), self.n)
    }

    /// Compute the restricted algebraic immunity of the Boolean function.
    ///
    /// # Arguments
    ///
    /// * `subset` - restriction set `S ⊆ {0,1}^n`
    ///
    /// # Returns
    ///
    /// Restricted Algebraic immunity of the Boolean function on S: `AI_S(f)`.
    ///
    /// # Example
    /// ```
    /// use algebraic_immunity::boolean_function::BooleanFunction;
    /// let t_table = vec![0,1,1,0];
    /// let bf = BooleanFunction::new(t_table);
    /// assert_eq!(bf.restricted_algebraic_immunity(vec![0,1]), 1);
    /// ```
    pub fn restricted_algebraic_immunity(&self, subset: Vec<usize>) -> usize {
        RestrictedAlgebraicImmunity::algebraic_immunity(self.t_table.clone(), subset, self.n)
    }
}
