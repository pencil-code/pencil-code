#
# USAGE:
#   make - make the full documentation using the Docker container
#   make fast - make the rst documentation (but no Python or Fortran autodocs) using the Docker container
#
default:
	docker compose up --build

fast:
	docker compose -f docker-compose-fast.yml up --build

