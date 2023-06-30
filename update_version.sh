#!/bin/bash

old_1=0
old_2=1
old_3=6

new_1=0
new_2=1
new_3=7

sed -i "s/${old_1}\.${old_2}\.${old_3}/${new_1}\.${new_2}\.${new_3}/g" $(find -not -path "./.git/*" -not -path "./target/*" -not -path "./Cargo.lock" -type f)