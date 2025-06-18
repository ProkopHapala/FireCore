from __future__ import annotations
import random
import math
import copy
import numpy as np

class Optimizer:
    """
    A generic Monte Carlo optimizer.
    """
    def __init__(self,
                 initial_solution: any,
                 evaluate_fitness: callable, # Args: (solution) -> tuple[float, str] (fitness, details_string)
                 mutation_callbacks: list[callable],
                 mutation_cumulative_probs: np.ndarray,
                 max_iterations: int = 1000,
                 pre_eval_range_checker: callable = None, # Args: (solution) -> tuple[bool, dict]
                 temperature_initial: float = 1.0,
                 temperature_decay: float = 0.995,
                 verbose: bool = True):
        """
        Initializes the optimizer.

        Parameters:
        -----------
        initial_solution : any
            The starting point for the optimization.
        evaluate_fitness : callable
            A function that takes a solution and returns a tuple: (fitness_score, details_string),
            where fitness_score is float (lower is better).
        mutation_callbacks : list[callable]
            A list of functions, where each function takes a solution and returns a new, mutated solution.
        mutation_cumulative_probs : np.ndarray
            A 1D NumPy array of cumulative probabilities corresponding to mutation_callbacks.
        pre_eval_range_checker : callable, optional # Args: (solution) -> tuple[bool, str]
            A function that takes a solution and returns (is_valid, details_string) for pre-evaluation checks.
        max_iterations : int
            The maximum number of iterations to run.
        temperature_initial : float
            Initial temperature for simulated annealing-like acceptance.
        temperature_decay : float
            Factor by which temperature decays each iteration.
        verbose : bool
            Whether to print progress.
        """
        self.current_solution = copy.deepcopy(initial_solution)
        self.evaluate_fitness = evaluate_fitness
        self.mutation_callbacks = mutation_callbacks
        self.mutation_cumulative_probs = mutation_cumulative_probs
        self.max_iterations = max_iterations
        self.pre_eval_range_checker = pre_eval_range_checker
        self.temperature = temperature_initial
        self.temperature_decay = temperature_decay
        self.verbose = verbose

        self.best_solution = copy.deepcopy(self.current_solution)
        self.current_fitness, self.current_details_str = self.evaluate_fitness(self.current_solution)
        self.best_fitness = self.current_fitness
        self.best_details_str = self.current_details_str
        
        self.history = [] # To store (iteration, fitness, details_string, temperature)

    def _apply_mutation(self, solution: any) -> any:
        """Selects and applies a mutation based on cumulative probabilities."""
        if not self.mutation_callbacks or self.mutation_cumulative_probs.size == 0:
            return copy.deepcopy(solution) # No mutations defined

        rand_val = random.random()
        mutation_idx = np.searchsorted(self.mutation_cumulative_probs, rand_val)
        
        if mutation_idx < len(self.mutation_callbacks):
            return self.mutation_callbacks[mutation_idx](solution)
        else: # Should not happen if probabilities are correctly normalized and sum to 1
            return copy.deepcopy(solution)

    def run(self):
        """
        Runs the optimization loop.
        """
        for i in range(self.max_iterations):
            new_solution_candidate = self._apply_mutation(self.current_solution)
            
            new_fitness = float('inf')
            new_details_str = "Rejected (pre-eval)"

            passed_pre_eval = True
            if self.pre_eval_range_checker:
                is_valid_pre, pre_eval_details_str = self.pre_eval_range_checker(new_solution_candidate)
                new_details_str = pre_eval_details_str # Store details string even if rejected
                if not is_valid_pre:
                    passed_pre_eval = False
            
            if passed_pre_eval:
                new_fitness, new_details_str = self.evaluate_fitness(new_solution_candidate)

            if new_fitness < self.current_fitness:
                self.current_solution = new_solution_candidate
                self.current_fitness = new_fitness
                self.current_details_str = new_details_str
                if new_fitness < self.best_fitness:
                    self.best_solution = copy.deepcopy(new_solution_candidate)
                    self.best_fitness = new_fitness
                    self.best_details_str = new_details_str
            elif self.temperature > 1e-6 and \
                 math.exp((self.current_fitness - new_fitness) / self.temperature) > random.random():
                self.current_solution = new_solution_candidate
                self.current_fitness = new_fitness
                self.current_details_str = new_details_str

            self.temperature *= self.temperature_decay
            self.history.append((i, self.current_fitness, self.current_details_str, self.temperature))

            if self.verbose and (i % (self.max_iterations // 20) == 0 or i == self.max_iterations -1) :
                print(f"Iter {i:5d}: BestFit={self.best_fitness:.3e} ({self.best_details_str}), CurrFit={self.current_fitness:.3e} ({self.current_details_str}), Temp={self.temperature:.2e}")
        
        if self.verbose:
            print(f"\nOptimization finished. Best fitness: {self.best_fitness:.4e} ({self.best_details_str})")
            print(f"Best basis structure: {self.best_solution}")
        return self.best_solution, self.best_fitness, self.history