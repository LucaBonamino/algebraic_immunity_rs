use itertools::Itertools;
use rayon::prelude::*;
use crate::vandermonde::{VanderMonde, str_ops, verify}; 
use lin_algebra::matrix::MatrixTrait;
use std::cmp::min;
use std::collections::HashSet;

#[derive(Clone)]
pub struct AlgebraicImmunity {
    truth_table: Vec<u8>
}

#[derive(Clone)]
pub struct RestrictedALgebraicImmunity {
    truth_table: Vec<u8>
}

pub trait AlgebraicImmunityTrait {
    // fn algebraic_immunity(truth_table: Vec<u8>, n: usize) -> usize;
    //fn find_min_annihilator(z: Vec<String>, e: Vec<String>, n: usize) -> Option<usize>;  
    fn generate_combinations(n: usize, r: usize) -> Vec<String> {
        let mut all_combinations = Vec::new();

        for k in 0..=r {
            for ones_positions in (0..n).combinations(k) {
                let mut binary_string = vec!['0'; n];
                for &pos in &ones_positions {
                    binary_string[pos] = '1';
                }
                let combination: String = binary_string.iter().rev().collect();
                all_combinations.push(combination);
            }
        }

        all_combinations
    } 
}


impl AlgebraicImmunity {

    pub fn new(truth_table: Vec<u8>) -> Self {

        let len = truth_table.len();
        assert!(len.is_power_of_two(), "Truth table length must be a power of two.");
        AlgebraicImmunity { truth_table }
    }

    fn compute_z(&self, n: usize) -> (Vec<String>, Vec<String>) {
        let mut true_idxs = Vec::new();
        let mut false_idxs = Vec::new();
     
        for i in 0..self.truth_table.len() {
           
            let bin_str = format!("{:0width$b}", i, width = n);
            if self.truth_table[i] == 1 {
                true_idxs.push(bin_str.clone());
            } else {
                false_idxs.push(bin_str.clone());
            }
            
        }

        (true_idxs, false_idxs)
    }

     /// Computes the algebraic immunity of a Boolean function of 'n' variables.
    /// 
    /// # Arguments
    /// 
    /// * 'truth_table' - A vecors of '1's and '0's representing the truth table of the Boolean function .
    /// * 'n' - Number of variables of the Boolean function.
    /// 
    /// # Returns
    /// 
    /// The algebraic immunity of the Boolean function as an integer.
    /// 
    /// # Examples
    /// 
    /// Algebraic immunity of constant function (1,1,1,1) -> the function f+1 (with truth table [0,0,0,0]) gets annihilates by g(x) = 1.
    /// ```
    /// use algebraic_immunity::ai::{AlgebraicImmunity, AlgebraicImmunityTrait};
    /// 
    /// let truth_table = vec![1,1,1,1];
    /// let n = 2;
    /// let ai = AlgebraicImmunity::algebraic_immunity(truth_table, n);
    /// assert_eq!(ai, 0);
    /// ```
    /// Functoin with algebraic immunity equal to 1.
    /// ```
    /// use algebraic_immunity::ai::{AlgebraicImmunity, AlgebraicImmunityTrait};
    /// 
    /// let truth_table = vec![0,1,0,0];
    /// let n = 2;
    /// let ai = AlgebraicImmunity::algebraic_immunity(truth_table, n);
    /// assert_eq!(ai, 1);
    /// ```
    pub fn algebraic_immunity(truth_table: Vec<u8>, n: usize) -> usize {
        let restricted_ai = Self::new(truth_table);
        let (z, z_c) = restricted_ai.compute_z(n);

        if z.is_empty() || z_c.is_empty() {
            return 0;
        }

        let r = (n+1) / 2;
        let e = Self::generate_combinations(n, r);

        let args = vec![
            (z.clone(), e.clone(), n),
            (z_c.clone(), e.clone(), n),
        ];

        let results: Vec<Option<usize>> = args
            .par_iter()
            .map(|(z, e, n)| {
                Self::find_min_annihilator(z.clone(), e.clone(), *n)
            })
            .collect();

        match results.into_iter().flatten().min() {
            Some(min_val) => min_val,
            None => 0,
        }
    }

    fn find_min_annihilator(
        mut z: Vec<String>,
        e: Vec<String>,
        n: usize
    ) -> Option<usize> {

        let max_number_of_monimials = e.len()-1;
        if max_number_of_monimials == 0{
            return None;
        }
        let size_support = z.len();
        if size_support < n+1{
            // If the cardinality of the support is smaller than D_1^n, the an annihiliator of degree d <= 1 must exist. if d was 0, 
            // it would hav been detcted ba the caller of this function. Therefore d = 1.
            return Some(1);
        }

        let mut vander_monde = VanderMonde::new(vec![
            vec![str_ops(&z[0], &e[0])]
        ]);
        

        let mut idx = 0;
        let mut i = 1;
        let mut operations: Vec<(usize, usize)> = vec![];

        let n_iters = min(size_support, max_number_of_monimials);

        while i < n_iters {

            vander_monde = vander_monde.compute_next(e[..=i].to_vec(), z[..=i].to_vec(), i, operations.clone());
            let (new_matrix, operations_i) = vander_monde.echelon_form();
            vander_monde = VanderMonde::from(new_matrix);

            if vander_monde.rank() < i + 1 {
                let kernel = vander_monde.kernel();
                // The kernel basis only contains maximum one element because of the algorithm design.
                let k = &kernel[0];

                let (vanish_on_z, vanish_index_opt) = verify(z[i + 1..].to_vec(), k.clone(), e[..=i].to_vec());
                if vanish_on_z {
                    return Some(hamming_weight(&e[i]));
                } else if let Some(vanish_index) = vanish_index_opt {
                    let new_index = i + vanish_index.0 + 1;
                    if new_index < z.len() {
                        z.swap(i + 1, new_index);
                    }
                }
            }

            i += 1;
            idx += 1;
            operations.extend(operations_i);
        }

        if (n_iters == size_support && size_support == max_number_of_monimials) || n_iters == max_number_of_monimials {
            // If the maximum number of iterations are reached, the algebraic immunity is ceil(n/2) - the hamming weight of the last monomial.
            if let Some(last) = e.last() {
                return Some(hamming_weight(&last));
            } else {
                return None;
            }
        } else if n_iters == size_support{
            // If all the elements of the support have been considered, at the next itaration, the matrix will not be squared anymore, and hence the rank of V_{n_iters+1} is not full anymore.
            return Some(hamming_weight(&e[idx+1]));
        }

        None

    }
 

}

impl AlgebraicImmunityTrait for AlgebraicImmunity {
    
}

impl RestrictedALgebraicImmunity {

    pub fn new(truth_table: Vec<u8>) -> Self {
        Self { truth_table }
    }

    fn compute_z(&self, subset: Vec<usize>, n: usize) -> (Vec<String>, Vec<String>, Vec<String>) {
        let mut true_idxs = Vec::new();
        let mut false_idxs = Vec::new();
        let mut s_bin = Vec::new();
        let s: HashSet<_> = subset.into_iter().collect();

        for i in 0..self.truth_table.len() {
            if !s.contains(&i) {
                continue;
            }

            let bin_str = format!("{:0width$b}", i, width = n);
            if self.truth_table[i] == 1 {
                true_idxs.push(bin_str.clone());
            } else {
                false_idxs.push(bin_str.clone());
            }
            s_bin.push(bin_str);
        }

        (true_idxs, false_idxs, s_bin)
    } 

    pub fn algebraic_immunity(truth_table: Vec<u8>, subset: Vec<usize>, n: usize) -> usize {
        let restricted_ai = Self::new(truth_table);
        let (z, z_c, s_bin) = restricted_ai.compute_z(subset, n);

        if z.is_empty() || z_c.is_empty() {
            return 0;
        }

        let e = Self::generate_combinations(n, n);

        let args = vec![
            (z.clone(), z_c.clone(), e.clone(), s_bin.clone()),
            (z_c.clone(), z.clone(), e.clone(), s_bin.clone()),
        ];

        let results: Vec<Option<usize>> = args
            .par_iter()
            .map(|(z, z_c, e, s_bin)| {
                RestrictedALgebraicImmunity::find_min_annihilator(z.clone(), z_c.clone(), e.clone(), s_bin.clone())
            })
            .collect();

        match results.into_iter().flatten().min() {
            Some(min_val) => min_val,
            None => 0,
        }
    }

    fn find_min_annihilator(
        mut z: Vec<String>,
        z_c: Vec<String>,
        mut e: Vec<String>,
        s: Vec<String>,
    ) -> Option<usize> {
        let mut vander_monde = VanderMonde::new(vec![
            vec![str_ops(&z[0], &e[0])]
        ]);

        let mut idx = 0;
        let mut i = 1;
        let mut operations: Vec<(usize, usize)> = vec![];

        let n_iters = z.len();

        while i < n_iters {
            let vander_monde_old = vander_monde.clone();

            vander_monde = vander_monde.compute_next(e[..=i].to_vec(), z[..=i].to_vec(), i, operations.clone());
            let (new_matrix, operations_i) = vander_monde.echelon_form();
            vander_monde = VanderMonde::from(new_matrix);

            if vander_monde.rank() < i + 1 {
                let kernel = vander_monde.kernel();
                let k = &kernel[0];

                let (vanish_on_z, vanish_index_opt) = verify(z[i + 1..].to_vec(), k.clone(), e[..=i].to_vec());
                if vanish_on_z {
                    let (vanish_on_s, _) = verify(z_c.clone(), k.clone(), e[..=i].to_vec());
                    if !vanish_on_s {
                        return Some(e[i].chars().filter(|c| *c == '1').count());
                    } else {
                        vander_monde = vander_monde_old;
                        e.remove(i);
                        continue;
                    }
                } else if let Some(vanish_index) = vanish_index_opt {
                    let new_index = i + vanish_index.0 + 1;
                    if new_index < z.len() {
                        z.swap(i + 1, new_index);
                    }
                }
            }

            i += 1;
            idx += 1;
            operations.extend(operations_i);
        }

        
        let mut vander_monde_s = VanderMonde::compute_vandermonde(s[..=idx].to_vec(), e[..=idx].to_vec());
        vander_monde_s = vander_monde_s.fill_rows(s[idx + 1..].to_vec(), e[..=idx].to_vec());

        let (vander_monde_s_reduced, mut operations_s) = vander_monde_s.echelon_form();
        let mut r_s = vander_monde_s_reduced.rank();

        if vander_monde.rank() < r_s {
            return Some(e[idx].chars().filter(|c| *c == '1').count());
        }

        i = idx + 1;
        let s_len = s.len();
        let mut vander_monde_s = VanderMonde::from(vander_monde_s_reduced);

        while r_s <= (s_len + 1) / 2 {
            if i >= e.len() {
                break;
            }

            vander_monde = vander_monde.construct_and_add_column(
                z.clone(),
                e[i].clone(),
                operations.clone()
            );

            vander_monde_s = vander_monde_s.construct_and_add_column(
                s.clone(),
                e[i].clone(),
                operations_s.clone()
            );

            let (vander_monde_s_new, ops_s) = vander_monde_s.echelon_form();
            vander_monde_s = VanderMonde::from(vander_monde_s_new);

            r_s = vander_monde_s.rank();

            if vander_monde.rank() < r_s {
                return Some(e[i].chars().filter(|c| *c == '1').count());
            }

            i += 1;
            operations_s.extend(ops_s);
        }


        None
    }
  
}

impl AlgebraicImmunityTrait for RestrictedALgebraicImmunity{}


fn hamming_weight(word: &str) -> usize{
    word.chars().filter(|c| *c == '1').count()
}


