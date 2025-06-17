from __future__ import annotations
import random
import math
import copy

class Optimizer:
    """
    A generic Monte Carlo optimizer.
    """
    def __init__(self,
                 initial_solution: any,
                 evaluate_fitness: callable, # Args: (solution) -> tuple[float, str] (fitness, details_string)
                 propose_mutation: callable, # Args: (solution) -> new_solution
                 max_iterations: int = 1000,
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
            A function that takes a solution and returns a tuple:
            (fitness_score, details_string), where fitness_score is float (lower is better).
        propose_mutation : callable
            A function that takes a solution and returns a new, mutated solution.
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
        self.propose_mutation = propose_mutation
        self.max_iterations = max_iterations
        self.temperature = temperature_initial
        self.temperature_decay = temperature_decay
        self.verbose = verbose

        self.best_solution = copy.deepcopy(self.current_solution)
        self.current_fitness, self.current_details_str = self.evaluate_fitness(self.current_solution)
        self.best_fitness = self.current_fitness
        self.best_details_str = self.current_details_str
        
        self.history = [] # To store (iteration, fitness, details_string, temperature)

    def run(self):
        """
        Runs the optimization loop.
        """
        for i in range(self.max_iterations):
            new_solution = self.propose_mutation(self.current_solution)
            new_fitness, new_details_str = self.evaluate_fitness(new_solution)

            if new_fitness < self.current_fitness:
                self.current_solution = new_solution
                self.current_fitness = new_fitness
                self.current_details_str = new_details_str
                if new_fitness < self.best_fitness:
                    self.best_solution = copy.deepcopy(new_solution)
                    self.best_fitness = new_fitness
                    self.best_details_str = new_details_str
            # Metropolis-like acceptance for worse solutions (simulated annealing)
            elif self.temperature > 1e-6 and \
                 math.exp((self.current_fitness - new_fitness) / self.temperature) > random.random():
                self.current_solution = new_solution
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