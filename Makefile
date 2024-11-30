.PHONY: init test

init:
	poetry install --no-dev

dev:
	poetry install
	poetry shell

test:
	python -m unittest test.test_basic
	python -m unittest test.test_show
	# python -m unittest test.test_align
	python -m unittest test.test_advance
