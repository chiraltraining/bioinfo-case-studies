# Quick Introduction to Python

## Printing Message

# To print something
print("Hello World!")

## Variables 
- What is Variables?
- Creating / Declaration of Variables 
- Naming Conventions 
- Restrictions / Reserved Words / Identifiers 

# Create a variable called `a` and store value
a = 10 
print(a)

# Create a variable called `b` and store value
b = 10.5 
print(b)

# Create a variable called `name` and store value
name = "Datastics Lab"
print(name)

## Reserved Keywords

import keyword
print(keyword.kwlist)

## Imports

# 'generic import' of math module
import math
math.sqrt(25)

# import a function
from math import sqrt
sqrt(25)    # no longer have to reference the module

# import multiple functions at once
from math import cos, floor

# import all functions in a module (generally discouraged)
from csv import *

# define an alias / nickname
import datetime as dt

# import numpy as np(nickname)
import numpy as np 

# show all functions in math module
print(dir(math))

## Data Types

**Determine the type of an object:**

type(2)

type(2.0)

type('two')

type(True)

type(None)

**Check if an object is of a given type:**

isinstance(2.0, int)

isinstance(2.0, (int, float))

**Convert an object to a given type:**

float(2)

int(2.9)

str(2.9)

**Zero, `None`, and empty containers are converted to `False`:**

bool(0)

bool(None)

bool('')    # empty string

bool([])    # empty list

bool({})    # empty dictionary

**Non-empty containers and non-zeros are converted to `True`:**

bool(2)

bool('two')

bool([2])

## Data Type Conversion

# float to integer  
a = 5.25 
int(a)

# integer to float  
b = 5
float(b) 

# integer to string  
c = 5
str(c)  

# string to integer  
d = '5'
str(d) 

# character to integer 
s = 'a' 
ord(s)

# integer to character 
i = 97 
chr(i)

## 3. Math

10 + 4

10 - 4

10 * 4

10 ** 4    # exponent

5 % 4      # modulo - computes the remainder

# Python 2: returns 2 (because both types are 'int')
# Python 3: returns 2.5
10 / 4

10 / float(4)

# force '/' in Python 2 to perform 'true division' (unnecessary in Python 3)
from __future__ import division

10 / 4     # true division

10 // 4    # floor division

## Comparisons and Boolean Operations

**Assignment statement:**

x = 5

**Comparisons:**

x > 3

x >= 3

x != 3

x == 5

**Boolean operations:**

5 > 3 and 6 > 3

5 > 3 or 5 < 3

not False

False or not False and True     # evaluation order: not, and, or

## Conditional Statements

# if statement
if x > 0:
    print('positive')

# if/else statement
if x > 0:
    print('positive')
else:
    print('zero or negative')

# if/elif/else statement
if x > 0:
    print('positive')
elif x == 0:
    print('zero')
else:
    print('negative')

# single-line if statement (sometimes discouraged)
if x > 0: print('positive')

# single-line if/else statement (sometimes discouraged), known as a 'ternary operator'
'positive' if x > 0 else 'zero or negative'

## Lists

- **List properties:** ordered, iterable, mutable, can contain multiple data types

# create an empty list (two ways)
empty_list = []
empty_list = list()

# create a list
simpsons = ['homer', 'marge', 'bart']

**Examine a list:**

# print element 0
simpsons[0]

len(simpsons)

**Modify a list (does not return the list):**

# append element to end
simpsons.append('lisa')
simpsons

# append multiple elements to end
simpsons.extend(['itchy', 'scratchy'])
simpsons

# insert element at index 0 (shifts everything right)
simpsons.insert(0, 'maggie')
simpsons

# search for first instance and remove it
simpsons.remove('bart')
simpsons

# remove element 0 and return it
simpsons.pop(0)

# remove element 0 (does not return it)
del simpsons[0]
simpsons

# replace element 0
simpsons[0] = 'krusty'
simpsons

# concatenate lists (slower than 'extend' method)
neighbors = simpsons + ['ned', 'rod', 'todd']
neighbors

**Find elements in a list:**

# counts the number of instances
simpsons.count('lisa')

# returns index of first instance
simpsons.index('itchy')

**List slicing:**

weekdays = ['mon', 'tues', 'wed', 'thurs', 'fri']

# element 0
weekdays[0]

# elements 0 (inclusive) to 3 (exclusive)
weekdays[0:3]

# starting point is implied to be 0
weekdays[:3]

# elements 3 (inclusive) through the end
weekdays[3:]

# last element
weekdays[-1]

# every 2nd element (step by 2)
weekdays[::2]

# backwards (step by -1)
weekdays[::-1]

# alternative method for returning the list backwards
list(reversed(weekdays))

**Sort a list in place (modifies but does not return the list):**

simpsons.sort()
simpsons

# sort in reverse
simpsons.sort(reverse=True)
simpsons

# sort by a key
simpsons.sort(key=len)
simpsons

**Return a sorted list (does not modify the original list):**

sorted(simpsons)

sorted(simpsons, reverse=True)

sorted(simpsons, key=len)

**Insert into an already sorted list, and keep it sorted:**

num = [10, 20, 40, 50]
from bisect import insort
insort(num, 30)
num

**Object references and copies:**

# create a second reference to the same list
same_num = num

# modifies both 'num' and 'same_num'
same_num[0] = 0
print(num)
print(same_num)

# copy a list (two ways)
new_num = num[:]
new_num = list(num)

**Examine objects:**

num is same_num    # checks whether they are the same object

num is new_num

num == same_num    # checks whether they have the same contents

num == new_num

## Tuples

- **Tuple properties:** ordered, iterable, immutable, can contain multiple data types
- Like lists, but they don't change size

# create a tuple directly
digits = (0, 1, 'two')

# create a tuple from a list
digits = tuple([0, 1, 'two'])

# trailing comma is required to indicate it's a tuple
zero = (0,)

**Examine a tuple:**

digits[2]

len(digits)

# counts the number of instances of that value
digits.count(0)

# returns the index of the first instance of that value
digits.index(1)

**Modify a tuple:**

# elements of a tuple cannot be modified (this would throw an error)
# digits[2] = 2

# concatenate tuples
digits = digits + (3, 4)
digits

**Other tuple operations:**

# create a single tuple with elements repeated (also works with lists)
(3, 4) * 2

# sort a list of tuples
tens = [(20, 60), (10, 40), (20, 30)]
sorted(tens)    # sorts by first element in tuple, then second element

# tuple unpacking
bart = ('male', 10, 'simpson')    # create a tuple
(sex, age, surname) = bart        # assign three values at once
print(sex)
print(age)
print(surname)

## Strings

- **String properties:** iterable, immutable

# convert another data type into a string
s = str(42)
s

# create a string directly
s = 'I like you'

**Examine a string:**

s[0]

len(s)

**String slicing is like list slicing:**

s[:6]

s[7:]

s[-1]

**Basic string methods (does not modify the original string):**

s.lower()

s.upper()

s.startswith('I')

s.endswith('you')

# checks whether every character in the string is a digit
s.isdigit()

# returns index of first occurrence, but doesn't support regex
s.find('like')

# returns -1 since not found
s.find('hate')

# replaces all instances of 'like' with 'love'
s.replace('like', 'love')

**Split a string:**

# split a string into a list of substrings separated by a delimiter
s.split(' ')

# equivalent (since space is the default delimiter)
s.split()

s2 = 'a, an, the'
s2.split(',')

**Join or concatenate strings:**

# join a list of strings into one string using a delimiter
stooges = ['larry', 'curly', 'moe']
' '.join(stooges)

# concatenate strings
s3 = 'The meaning of life is'
s4 = '42'
s3 + ' ' + s4

**Remove whitespace from the start and end of a string:**

s5 = '  ham and cheese  '
s5.strip()

**String substitutions:**

# old way
'raining %s and %s' % ('cats', 'dogs')

# new way
'raining {} and {}'.format('cats', 'dogs')

# new way (using named arguments)
'raining {arg1} and {arg2}'.format(arg1='cats', arg2='dogs')

**String formatting ([more examples](https://mkaz.tech/python-string-format.html)):**

# use 2 decimal places
'pi is {:.2f}'.format(3.14159)

**Normal strings versus raw strings:**

# normal strings allow for escaped characters
print('first line\nsecond line')

# raw strings treat backslashes as literal characters
print(r'first line\nfirst line')

## Dictionaries

- **Dictionary properties:** unordered, iterable, mutable, can contain multiple data types
- Made of key-value pairs
- Keys must be unique, and can be strings, numbers, or tuples
- Values can be any type

# create an empty dictionary (two ways)
empty_dict = {}
empty_dict = dict()

# create a dictionary (two ways)
family = {'dad':'homer', 'mom':'marge', 'size':6}
family = dict(dad='homer', mom='marge', size=6)
family

# convert a list of tuples into a dictionary
list_of_tuples = [('dad', 'homer'), ('mom', 'marge'), ('size', 6)]
family = dict(list_of_tuples)
family

**Examine a dictionary:**

# pass a key to return its value
family['dad']

# return the number of key-value pairs
len(family)

# check if key exists in dictionary
'mom' in family

# dictionary values are not checked
'marge' in family

# returns a list of keys (Python 2) or an iterable view (Python 3)
family.keys()

# returns a list of values (Python 2) or an iterable view (Python 3)
family.values()

# returns a list of key-value pairs (Python 2) or an iterable view (Python 3)
family.items()

**Modify a dictionary (does not return the dictionary):**

# add a new entry
family['cat'] = 'snowball'
family

# edit an existing entry
family['cat'] = 'snowball ii'
family

# delete an entry
del family['cat']
family

# dictionary value can be a list
family['kids'] = ['bart', 'lisa']
family

# remove an entry and return the value
family.pop('dad')

# add multiple entries
family.update({'baby':'maggie', 'grandpa':'abe'})
family

**Access values more safely with `get`:**

family['mom']

# equivalent to a dictionary lookup
family.get('mom')

# this would throw an error since the key does not exist
# family['grandma']

# return None if not found
family.get('grandma')

# provide a default return value if not found
family.get('grandma', 'not found')

**Access a list element within a dictionary:**

family['kids'][0]

family['kids'].remove('lisa')
family

**String substitution using a dictionary:**

'youngest child is %(baby)s' % family

## Sets

- **Set properties:** unordered, iterable, mutable, can contain multiple data types
- Made of unique elements (strings, numbers, or tuples)
- Like dictionaries, but with keys only (no values)

# create an empty set
empty_set = set()

# create a set directly
languages = {'python', 'r', 'java'}

# create a set from a list
snakes = set(['cobra', 'viper', 'python'])

**Examine a set:**

len(languages)

'python' in languages

**Set operations:**

# intersection
languages & snakes

# union
languages | snakes

# set difference
languages - snakes

# set difference
snakes - languages

**Modify a set (does not return the set):**

# add a new element
languages.add('sql')
languages

# try to add an existing element (ignored, no error)
languages.add('r')
languages

# remove an element
languages.remove('java')
languages

# try to remove a non-existing element (this would throw an error)
# languages.remove('c')

# remove an element if present, but ignored otherwise
languages.discard('c')
languages

# remove and return an arbitrary element
languages.pop()

# remove all elements
languages.clear()
languages

# add multiple elements (can also pass a set)
languages.update(['go', 'spark'])
languages

**Get a sorted list of unique elements from a list:**

sorted(set([9, 0, 2, 1, 0]))

## Defining Functions

**Define a function with no arguments and no return values:**

def print_text():
    print('this is text')

# call the function
print_text()

**Define a function with one argument and no return values:**

def print_this(x):
    print(x)

# call the function
print_this(3)

# prints 3, but doesn't assign 3 to n because the function has no return statement
n = print_this(3)

**Define a function with one argument and one return value:**

def square_this(x):
    return x**2

# include an optional docstring to describe the effect of a function
def square_this(x):
    """Return the square of a number."""
    return x**2

# call the function
square_this(3)

# assigns 9 to var, but does not print 9
var = square_this(3)

**Define a function with two 'positional arguments' (no default values) and one 'keyword argument' (has a default value):**


def calc(a, b, op='add'):
    if op == 'add':
        return a + b
    elif op == 'sub':
        return a - b
    else:
        print('valid operations are add and sub')

# call the function
calc(10, 4, op='add')

# unnamed arguments are inferred by position
calc(10, 4, 'add')

# default for 'op' is 'add'
calc(10, 4)

calc(10, 4, 'sub')

calc(10, 4, 'div')

**Use `pass` as a placeholder if you haven't written the function body:**

def stub():
    pass

**Return two values from a single function:**

def min_max(nums):
    return min(nums), max(nums)

# return values can be assigned to a single variable as a tuple
nums = [1, 2, 3]
min_max_num = min_max(nums)
min_max_num

# return values can be assigned into multiple variables using tuple unpacking
min_num, max_num = min_max(nums)
print(min_num)
print(max_num)

## Anonymous (Lambda) Functions

- Primarily used to temporarily define a function for use by another function

# define a function the "usual" way
def squared(x):
    return x**2

# define an identical function using lambda
squared = lambda x: x**2

**Sort a list of strings by the last letter:**

# without using lambda
simpsons = ['homer', 'marge', 'bart']
def last_letter(word):
    return word[-1]
sorted(simpsons, key=last_letter)

# using lambda
sorted(simpsons, key=lambda word: word[-1])

[<a href="#Python-Quick-Guide">Back to top</a>]

## For Loops and While Loops

**`range` returns a list of integers (Python 2) or a sequence (Python 3):**

# includes the start value but excludes the stop value
range(0, 3)

# default start value is 0
range(3)

# third argument is the step value
range(0, 5, 2)

# Python 2 only: use xrange to create a sequence rather than a list (saves memory)
xrange(100, 100000, 5)

**`for` loops:**

# not the recommended style
fruits = ['apple', 'banana', 'cherry']
for i in range(len(fruits)):
    print(fruits[i].upper())

# recommended style
for fruit in fruits:
    print(fruit.upper())

# iterate through two things at once (using tuple unpacking)
family = {'dad':'homer', 'mom':'marge', 'size':6}
for key, value in family.items():
    print(key, value)

# use enumerate if you need to access the index value within the loop
for index, fruit in enumerate(fruits):
    print(index, fruit)

**`for`/`else` loop:**

for fruit in fruits:
    if fruit == 'banana':
        print('Found the banana!')
        break    # exit the loop and skip the 'else' block
else:
    # this block executes ONLY if the for loop completes without hitting 'break'
    print("Can't find the banana")

**`while` loop:**

count = 0
while count < 5:
    print('This will print 5 times')
    count += 1    # equivalent to 'count = count + 1'

## Comprehensions

**List comprehension:**

# for loop to create a list of cubes
nums = [1, 2, 3, 4, 5]
cubes = []
for num in nums:
    cubes.append(num**3)
cubes

# equivalent list comprehension
cubes = [num**3 for num in nums]
cubes

# for loop to create a list of cubes of even numbers
cubes_of_even = []
for num in nums:
    if num % 2 == 0:
        cubes_of_even.append(num**3)
cubes_of_even

# equivalent list comprehension
# syntax: [expression for variable in iterable if condition]
cubes_of_even = [num**3 for num in nums if num % 2 == 0]
cubes_of_even

# for loop to cube even numbers and square odd numbers
cubes_and_squares = []
for num in nums:
    if num % 2 == 0:
        cubes_and_squares.append(num**3)
    else:
        cubes_and_squares.append(num**2)
cubes_and_squares

# equivalent list comprehension (using a ternary expression)
# syntax: [true_condition if condition else false_condition for variable in iterable]
cubes_and_squares = [num**3 if num % 2 == 0 else num**2 for num in nums]
cubes_and_squares

# for loop to flatten a 2d-matrix
matrix = [[1, 2], [3, 4]]
items = []
for row in matrix:
    for item in row:
        items.append(item)
items

# equivalent list comprehension
items = [item for row in matrix
              for item in row]
items

**Set comprehension:**

fruits = ['apple', 'banana', 'cherry']
unique_lengths = {len(fruit) for fruit in fruits}
unique_lengths

**Dictionary comprehension:**

fruit_lengths = {fruit:len(fruit) for fruit in fruits}
fruit_lengths

fruit_indices = {fruit:index for index, fruit in enumerate(fruits)}
fruit_indices

## Map and Filter

**`map` applies a function to every element of a sequence and returns a list (Python 2) or iterator (Python 3):**

simpsons = ['homer', 'marge', 'bart']
map(len, simpsons)

# equivalent list comprehension
[len(word) for word in simpsons]

map(lambda word: word[-1], simpsons)

# equivalent list comprehension
[word[-1] for word in simpsons]

**`filter` returns a list (Python 2) or iterator (Python 3) containing the elements from a sequence for which a condition is `True`:**

nums = range(5)
filter(lambda x: x % 2 == 0, nums)

# equivalent list comprehension
[num for num in nums if num % 2 == 0]

## Randomness 

# import random module 
import random

# generate random numbers 
random.random() 

# generate random numbers 
random.random()

# generate random numbers from start, end 
random.randint(1, 20)

# generate random numbers from start, end 
random.randint(1, 20)

# generate random numbers using randrange(start, step, stop) 
random.randrange(1, 20, 2)

# generate random numbers using randrange(start, step, stop) 
random.randrange(1, 20, 2)

# randomly choice: random.choice([sequence]) 
random.choice([1, 2, 3, 4, 5, 6])

# randomly choice: random.choice([sequence]) 
random.choice([1, 2, 3, 4, 5, 6])

# randomly choice: random.choice([sequence]) 
random.choice(["H", "T"])

# randomly choice: random.choice([sequence]) 
random.choice(["H", "T"])

## Working with Files 

# store file in filename variable 
filename = "./data/input.txt"

# open file 
open(filename)

# read file content using for loop 
for line in open(filename): 
    print(line)

# removed special character: \n and \t
for line in open(filename): 
    line = line.rstrip() 
    print(line)

# split: returns a list of each line 
for line in open(filename):
    line = line.rstrip().split(" ")
    print(line)
    

# write a file 
F = open("./data/output.txt", "w") # w for writing mood 
F.write("Python\n")

# close file 
F.close() 

# read a file using with statement: open and close  
with open(filename, "r") as F: # r means reading mood 
    for line in F: 
        line = line.rstrip() 
        print(line)

# write a custom function for reading file 
def readFile(inputfile): 
    with open(inputfile, "r") as F: 
        for line in F: 
            line = line.rstrip()
            print(line)

# read file using readFile 
filename = "./data/input.txt"
readFile(filename)

## Common Mistakes and Errors

# Create a list 
L = [2,4,6]

# common error-1: IndexError
L[4]

# check length 
len(L)

# common error-2: AttributeError
L.add(8)

# solution 
L.append(8)
L 

# Create a dictionary 
D = {1: "one", 2:"two"}

# keys 
D.keys() 

# common error-3: KeyError
D[0]

# common error-4: TypeError 
"strings" + 9

# common error-5: IndentationError
def rsum(n): 
    rsum = 0 
    for k in range(n): 
        rsum += k 
        return rsum
rsum(12)

## References and Further Reading
- http://www.dataschool.io/
- https://www.python.org/doc/essays/blurb/
- https://dev.to/duomly/10-reasons-why-learning-python-is-still-a-great-idea-5abh
- https://www.stat.washington.edu/~hoytak/blog/whypython.html
- https://www.programiz.com/python-programming