.PHONY: init test

init:
	pipenv install

test:
	python -m unittest test.test_basic
	# python -m unittest test.test_show
	# python -m unittest test.test_align
	python -m unittest test.test_advance
