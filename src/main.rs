#![feature(isqrt)]
use primal::Sieve;
use std::collections::HashSet;

fn sig(n: usize) -> u16 {
    let mut out = 0;
    let mut work = n;
    while work > 0 {
        let digit = work % 10;
        out |= 1 << digit;
        work /= 10;
    }
    out
}

fn melissa(limit: usize) {
    let sieve = Sieve::new(limit + 1000);
    let mut primes = sieve.primes_from(2);
    let mut max_streak = 0;
    let mut streak = 0;
    let mut current_prime = primes.next().expect("More primes");
    'main: for n in 2..=limit {
        if current_prime < n {
            current_prime = primes.next().expect("More primes");
        }
        assert!(current_prime >= n);
        if current_prime - n + streak <= max_streak {
            streak = 0;
            continue;
        }
        let next_five = ((n + 4) / 10) * 10 + 5;
        assert!(next_five >= n);
        if next_five - n + streak <= max_streak {
            streak = 0;
            continue;
        }
        let n_sig = sig(n);
        if n_sig.count_ones() >= 5 {
            streak = 0;
            continue;
        }
        let factors = sieve.factor(n).expect("sqrt can't break");
        let mut prod = factors.iter().fold(1, |acc, f| acc * (f.1 + 1));
        // Never need most of the largest prime in a single divisor.
        // Can leave that portion in work.
        let lk = factors.last().expect("Final").1;
        prod /= lk + 1;
        prod *= lk / 2 + 1;
        // Targets are divisors of n that overlap with n
        let mut targets = HashSet::new();
        targets.insert(n);
        for mut exps in 1..prod {
            let mut factor = 1;
            for (p, k) in factors.iter() {
                let exp = exps % (k + 1);
                factor *= p.pow(exp as u32);
                exps /= k + 1;
            }
            let i_sig = sig(factor);
            if i_sig & n_sig == 0 {
                let mut staging = vec![];
                for &target in &targets {
                    if target < factor {
                        continue;
                    }
                    let mut work = target;
                    while work % factor == 0 {
                        work /= factor;
                        if n_sig & sig(work) == 0 {
                            streak += 1;
                            let low = 6777666;
                            let run = 7;
                            if low <= n && n < low + run {
                                println!(
                                    "{n} {factor} {work} {} {}",
                                    factor * work,
                                    n / (factor * work)
                                );
                            }
                            if streak > max_streak {
                                max_streak = streak;
                                println!("{streak}: {}-{n}", n - streak + 1);
                            }
                            continue 'main;
                        }
                        staging.push(work);
                    }
                }
                targets.extend(staging);
            }
        }
        streak = 0;
    }
}

fn old_melissa(limit: usize) {
    let mut table: Vec<HashSet<u16>> = vec![HashSet::new()];
    let mut last = 0;
    let mut streak = 0;
    let mut max_streak = 0;
    let sieve = Sieve::new(limit.isqrt());
    for n in 1..=limit {
        let mut n_reps = HashSet::new();
        let n_sig = sig(n);
        n_reps.insert(n_sig);
        //let factors = Factorization::run(n).prime_factor_repr();
        let factors = sieve.factor(n).expect("sqrt can't break");
        let mut prod = factors.iter().fold(1, |acc, f| acc * (f.1 + 1));
        if let Some((_lp, lk)) = factors.last() {
            prod /= lk + 1;
            prod *= lk / 2 + 1;
        }
        for mut exps in 0..prod {
            let mut factor = 1;
            for (p, k) in factors.iter() {
                let exp = exps % (k + 1);
                factor *= p.pow(exp as u32);
                exps /= k + 1;
            }
            let oth_factor = n / factor;
            if factor > 1 && factor <= oth_factor {
                let i_reps = &table[factor as usize];
                let j_reps = &table[oth_factor as usize];
                for i_rep in i_reps {
                    for j_rep in j_reps {
                        let new_rep = i_rep | j_rep;
                        if new_rep & n_sig == 0 && n != last {
                            if n == last + 1 {
                                streak += 1;
                                if streak > max_streak {
                                    max_streak = streak;
                                    println!("{streak}: {}-{n}", n - streak + 1);
                                }
                            } else {
                                streak = 1;
                            }
                            last = n;
                        }
                        n_reps.insert(i_rep | j_rep);
                    }
                }
            }
        }
        table.push(n_reps);
    }
}
fn main() {
    let alg = std::env::args()
        .nth(1)
        .expect("alg given")
        .parse()
        .expect("alg parsed");
    let limit = std::env::args()
        .nth(2)
        .expect("limit given")
        .parse()
        .expect("limit parsed");
    match alg {
        0 => melissa(limit),
        1 => old_melissa(limit),
        _ => unimplemented!(),
    }
}
