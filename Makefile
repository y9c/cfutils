.PHONY: init test

init:
	peotry install --no-dev

dev:
	peotry install
	peotry shell

test:
	python -m unittest test.test_basic
	python -m unittest test.test_show
	# python -m unittest test.test_align
	python -m unittest test.test_advance
