#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 09:45:28 2019

@author: rubens
"""

"""
faça um programa que leia um número qualquer
e mostre o seu fatorial"""


from math import factorial
f=factorial(5)
print (f)

""" outra forma """

fat=1
n=5
#n= int(input('digite um número para calcular o fatorial: '))
c=n
while c>0:
    print('{}'.format(c), end='')
    print(' x ' if c > 1 else ' = ',end='')
    fat=fat*c
    c=c-1
print('{}'.format(f), end='')



"""com for

fat=1
n=5
for i in range (1,n+1):
    fat=fat*i
print (fat)

"""