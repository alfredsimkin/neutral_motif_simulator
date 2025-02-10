'''
trying to learn to use python's built in filecmp module to compare the identity of two files
'''

import filecmp

#file1=open('file1')
#file2=open('file2')

print filecmp.cmp('file1', 'file2')
