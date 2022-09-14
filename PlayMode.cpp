#include "PlayMode.hpp"

#include "LitColorTextureProgram.hpp"

#include "DrawLines.hpp"
#include "Mesh.hpp"
#include "Load.hpp"
#include "gl_errors.hpp"
#include "data_path.hpp"

//#include "../nest-libs/windows/glm/include/glm/gtc/type_ptr.hpp"
#include <glm/gtc/type_ptr.hpp>
//#include "../nest-libs/windows/glm/include/glm/glm.hpp"
#include <glm/glm.hpp>
//#include "../nest-libs/windows/glm/include/glm/gtc/quaternion.hpp"
#include <glm/gtc/quaternion.hpp>

#include <random>

//#include "../nest-libs/windows/SDL2/include/SDL.h"
#include <SDL.h>
//#include "../nest-libs/windows/glm/include/glm/gtx/color_space.hpp"
#include <glm/gtx/color_space.hpp>

GLuint cannon_meshes_for_lit_color_texture_program = 0;
Load< MeshBuffer > cannon_meshes(LoadTagDefault, []() -> MeshBuffer const * {
	MeshBuffer const *ret = new MeshBuffer(data_path("cannon.pnct"));
	cannon_meshes_for_lit_color_texture_program = ret->make_vao_for_program(lit_color_texture_program->program);
	return ret;
});

Load< Scene > cannon_scene(LoadTagDefault, []() -> Scene const * {
	return new Scene(data_path("cannon.scene"), [&](Scene &scene, Scene::Transform *transform, std::string const &mesh_name){
		// Cannon parts are kept in place; all others will need to be instantiated dynamically
		if (transform->name == "Pole" || transform->name == "Holder" || transform->name == "Cannon" || transform->name == "Nozzle") {
			Mesh const& mesh = cannon_meshes->lookup(mesh_name);

			scene.drawables.emplace_back(transform);
			Scene::Drawable& drawable = scene.drawables.back();

			drawable.pipeline = lit_color_texture_program_pipeline;

			drawable.pipeline.vao = cannon_meshes_for_lit_color_texture_program;
			drawable.pipeline.type = mesh.type;
			drawable.pipeline.start = mesh.start;
			drawable.pipeline.count = mesh.count;
		}
	});
});

glm::ivec2 PlayMode::getTilePos(GroundTile* tile) {
	size_t index = tile - map.data();
	int x = (int)index % map_diameter - map_radius;
	int y = (int)index / map_diameter - map_radius;
	return glm::ivec2(x, y);
}

glm::vec2 PlayMode::getTileCoords(GroundTile* tile) {
	return (glm::vec2)getTilePos(tile) * tile_size;
}

PlayMode::GroundTile* PlayMode::getTileAtPos(glm::ivec2 pos) {
	int x = pos[0];
	int y = pos[1];

	if (!isInRing(glm::ivec2(x, y))) {
		return nullptr;
	}

	int index = (y + map_radius) * map_diameter + (x + map_radius);
	if (map[index].exists && !map[index].falling) {
		return &map[index];
	}
	else {
		return nullptr;
	}
}

PlayMode::GroundTile* PlayMode::getTileAtCoords(glm::vec2 coords) {
	int x = (int)round(coords[0] / tile_size);
	int y = (int)round(coords[1] / tile_size);
	return getTileAtPos(glm::ivec2(x, y));
}

bool PlayMode::isInRing(glm::ivec2 pos) {
	float r = (float)sqrt(pos[0] * pos[0] + pos[1] * pos[1]);
	return r < map_radius && r >= hole_radius;
}

bool PlayMode::isInRing(GroundTile* tile) {
	return isInRing(getTilePos(tile));
}

void PlayMode::destroyTile(GroundTile* tile) {
	// Tile disappears
	tile->exists = false;
	tile->transform->position[2] = -100.f;

	// Incur score penalty for meeple death
	if (tile->meeple) {
		score += score_kill;
		perfect_score += score_healthy;
	}

	// Get tile grid position
	glm::ivec2 pos = getTilePos(tile);

	// All adjacent or diagonal tiles, in clockwise order from top left
	std::array<glm::ivec2, 8> loop = {
		glm::ivec2(-1, 1),
		glm::ivec2(0, 1),
		glm::ivec2(1, 1),
		glm::ivec2(1, 0),
		glm::ivec2(1, -1),
		glm::ivec2(0, -1),
		glm::ivec2(-1, -1),
		glm::ivec2(-1, 0)
	};

	// Partition neighbors into "chunks" of contiguous neighbors, and label them accordingly
	uint8_t chunk = 1;
	std::vector<GroundTile*> chunk_roots;
	for (size_t i = 0; i < loop.size(); i++) {
		GroundTile* other = getTileAtPos(pos + loop[i]);
		if (other != nullptr) {
			other->chunk = chunk;
			if (chunk_roots.empty() || chunk_roots.back()->chunk != chunk) {
				chunk_roots.push_back(other);
			}

			// If this is the last neighbor and it's right next to the first, they should be in the same chunk
			if (i == loop.size() - 1 && getTilePos(chunk_roots[0]) == pos + loop[0]) {
				chunk_roots.erase(chunk_roots.begin());
			}
		}
		else {
			chunk++;
		}
	}


	// If destroying this block didn't partition its neighbors, nothing needs to be done
	if (chunk_roots.size() <= 1) {
		return;
	}

	// Find the biggest chunk (by population)
	size_t max_population = 0;
	uint8_t biggest_chunk = 0;
	for (size_t i = 0; i < chunk_roots.size(); i++) {
		// This chunk was part of an earlier chunk
		if (chunk_roots[i]->visited_for == tile) {
			continue;
		}

		// Breadth-first search from chunk root, counting tiles visited and healthy/sick meeples
		size_t healthy = 0;
		size_t sick = 0;
		std::deque<GroundTile*> to_visit;
		chunk_roots[i]->visited_for = tile;
		to_visit.push_back(chunk_roots[i]);
		while (!to_visit.empty()) {
			GroundTile* visiting = to_visit.front();
			to_visit.pop_front();
			glm::ivec2 vpos = getTilePos(visiting);
			for (size_t j = 0; j < directions.size(); j++) {
				GroundTile* neighbor = getTileAtPos(vpos + directions[j]);
				if (neighbor != nullptr && neighbor->visited_for != tile) {
					neighbor->chunk = chunk_roots[i]->chunk;
					neighbor->visited_for = tile;
					to_visit.push_back(neighbor);

					if (neighbor->meeple) {
						if (neighbor->meeple->sick) {
							sick++;
						}
						else {
							healthy++;
						}
					}
				}
			}
		}

		size_t population = sick + healthy;

		// Update biggest chunk if applicable
		// Also, even if it is biggest, always drop "pure" tiles of all sick or all healthy people
		if (population > max_population && sick > 0 && healthy > 0) {
			max_population = population;
			biggest_chunk = chunk_roots[i]->chunk;
		}
	}

	// Destroy all smaller chunks
	for (size_t i = 0; i < chunk_roots.size(); i++) {
		if (chunk_roots[i]->chunk != biggest_chunk) {
			dropChunk(chunk_roots[i]);
		}
	}
}

void PlayMode::dropChunk(GroundTile* root) {
	size_t healthy = 0;
	size_t sick = 0;
	std::deque<GroundTile*> to_visit;
	root->falling = true;
	to_visit.push_back(root);
	while (!to_visit.empty()) {
		GroundTile* visiting = to_visit.front();
		to_visit.pop_front();
		glm::ivec2 vpos = getTilePos(visiting);
		for (size_t j = 0; j < directions.size(); j++) {
			GroundTile* neighbor = getTileAtPos(vpos + directions[j]);
			if (neighbor != nullptr) {
				neighbor->falling = true;
				to_visit.push_back(neighbor);

				if (neighbor->meeple) {
					if (neighbor->meeple->sick) {
						sick++;
					}
					else {
						healthy++;
					}
				}
			}
		}
	}

	// Calculate score from this chunk
	if (sick == 0) {
		score += (int)healthy * score_healthy;
	}
	else {
		score += (int)(healthy + sick) * score_sick;
	}
	perfect_score += (int)healthy * score_healthy + (int)sick * score_sick;
}

void PlayMode::setMesh(Scene::Drawable* drawable, std::string mesh_name) {
	// Assign the mesh with the given name to the given drawable
	Mesh const& mesh = cannon_meshes->lookup(mesh_name);
	drawable->pipeline = lit_color_texture_program_pipeline;
	drawable->pipeline.vao = cannon_meshes_for_lit_color_texture_program;
	drawable->pipeline.type = mesh.type;
	drawable->pipeline.start = mesh.start;
	drawable->pipeline.count = mesh.count;
}

void PlayMode::infect(Meeple* meeple, float exposure) {
	float original_infection = meeple->infection;
	
	// Nothing to do if they're already sick
	if (original_infection >= infection_threshold) {
		return;
	}

	// Increase infection counter
	meeple->infection += exposure;

	// If infection crossed the final threshold, meeple is sick
	if (meeple->infection >= infection_threshold && !meeple->sick) {
		meeple->sick = true;
		perfect_score += score_healthy - score_sick;
	}
	
	// Set meeple mesh based on level of infection
	for (size_t i = 0; i < num_thresholds; i++) {
		float threshold = infection_threshold * (1.f - (float)i / (float)num_thresholds);
		int threshold_num = (int)num_thresholds - (int)i;
		if (original_infection < threshold && meeple->infection >= threshold) {
			setMesh(meeple->drawable, "Meeple_" + std::to_string(threshold_num));
			break;
		}
	}
}

void PlayMode::reset(int level) {
	// Reset score
	score = 0;
	perfect_score = 0;

	// Reset target
	target = glm::vec3(0.f, -1.5f, 0.f);

	// Clear meeples
	for (size_t i = 0; i < meeples.size(); i++) {
		const Scene::Drawable drawable = *meeples[i]->drawable;
		for (auto it = scene.drawables.begin(); it != scene.drawables.end(); it++) {
			if (it->pipeline.start == drawable.pipeline.start && it->transform == drawable.transform) {
				scene.drawables.erase(it);
				break;
			}
		}
		delete drawable.transform;
		delete meeples[i];
	}
	meeples.clear();

	// Delete clouds
	while (!clouds.empty()) {
		deleteCloud(0);
	}

	// Initialize map
	for (int i = 0; i < map_diameter; i++) {
		for (int j = 0; j < map_diameter; j++) {
			int x = j - map_radius;
			int y = i - map_radius;
			size_t index = i * map_diameter + j;
			if (isInRing(glm::ivec2(x, y))) {
				// Add ground tile to map
				map[index].transform->position = glm::vec3(x * tile_size, y * tile_size, 0);
				map[index].exists = true;
				map[index].falling = false;
				map[index].meeple = nullptr;
				map[index].visited_for = nullptr;
				map[index].velocity = glm::vec3(0.f, 0.f, 0.f);
			}
			else {
				// No tile here
				map[index].exists = false;
			}
		}
	}

	// Set difficulty
	if (level > 0) {
		setLevel(level);
	}

	// Populate board with meeples
	placeMeeples();

	// Reset perfect score
	perfect_score = 0;
}

void PlayMode::setLevel(int level) {
	start_meeples = 50 + level * 10;
}

void PlayMode::placeMeeples() {
	for (size_t i = 0; i < start_meeples; i++) {
		size_t max_attempts = 100;
		for (size_t attempt = 0; attempt < max_attempts; attempt++) {
			size_t index = rand() % map.size();
			if (map[index].exists && !map[index].meeple) {
				//scene.transforms.emplace_back();
				Scene::Transform* transform = new Scene::Transform;
				transform->position = map[index].transform->position;
				transform->name = "Meeple_0";

				scene.drawables.emplace_back(transform);
				setMesh(&scene.drawables.back(), "Meeple_0");

				meeples.push_back(new Meeple);
				meeples.back()->drawable = &scene.drawables.back();
				meeples.back()->tile = &map[index];
				float theta = (float)rand() / (float)RAND_MAX * 2.f * float(M_PI);
				meeples.back()->dir = glm::vec2(cosf(theta), sinf(theta));
				meeples.back()->jump_wait_timer = (float)rand() / (float)RAND_MAX * jump_wait_duration;
				meeples.back()->sick = false;
				meeples.back()->infection = 0;
				meeples.back()->jumping_from = nullptr;
				meeples.back()->jumping_to = nullptr;
				map[index].meeple = meeples.back();

				break;
			}
		}
	}

	size_t num_sick = (size_t)((float)meeples.size() * start_sick);
	for (size_t i = 0; i < num_sick; i ++) {
		infect(meeples[i], infection_threshold);
	}
}

void PlayMode::deleteCloud(size_t i) {
	const Scene::Drawable drawable = *clouds[i]->drawable;
	for (auto it = scene.drawables.begin(); it != scene.drawables.end(); it++) {
		if (it->pipeline.start == drawable.pipeline.start && it->transform == drawable.transform) {
			scene.drawables.erase(it);
			break;
		}
	}
	delete drawable.transform;
	delete clouds[i];
	clouds.erase(clouds.begin() + i);
}

PlayMode::PlayMode() : scene(*cannon_scene) {
	// Set random seed for a different experience every time
	srand((int)time(nullptr));

	// Initialize map
	for (int i = 0; i < map_diameter; i++) {
		for (int j = 0; j < map_diameter; j++) {
			int x = j - map_radius;
			int y = i - map_radius;
			size_t index = i * map_diameter + j;
			if (isInRing(glm::ivec2(x, y))) {
				// Create ground tile transform
				scene.transforms.emplace_back();
				scene.transforms.back().position = glm::vec3(x * tile_size, y * tile_size, 0);
				scene.transforms.back().name = "Ground";

				// Create drawable for ground tile
				scene.drawables.emplace_back(&scene.transforms.back());
				if (index % 2 == 0) {
					setMesh(&scene.drawables.back(), "Ground_0");
				}
				else {
					setMesh(&scene.drawables.back(), "Ground_1");
				}
				
				// Add ground tile to map
				map[index].transform = &scene.transforms.back();
				map[index].exists = true;
			}
			else {
				// No tile here
				map[index].exists = false;
			}
		}
	}

	// Set level to determine meeple population, then create meeples
	setLevel(1);
	placeMeeples();

	// Reset perfect score (placeMeeples modifies it by infecting meeples)
	perfect_score = 0;

	// Get cannon part transforms
	for (auto& transform : scene.transforms) {
		if (transform.name == "Pole") pole = &transform;
		else if (transform.name == "Cannon") cannon = &transform;
		else if (transform.name == "Nozzle") nozzle = &transform;
	}

	//get pointer to camera for convenience:
	if (scene.cameras.size() != 1) throw std::runtime_error("Expecting scene to have exactly one camera, but it has " + std::to_string(scene.cameras.size()));
	camera = &scene.cameras.front();
}

PlayMode::~PlayMode() {
}

bool PlayMode::handle_event(SDL_Event const &evt, glm::uvec2 const &window_size) {

	if (evt.type == SDL_KEYDOWN) {
		if (evt.key.keysym.sym == SDLK_ESCAPE) {
			SDL_SetRelativeMouseMode(SDL_FALSE);
			return true;
		} else if (evt.key.keysym.sym == SDLK_a) {
			left.downs += 1;
			left.pressed = true;
			return true;
		} else if (evt.key.keysym.sym == SDLK_d) {
			right.downs += 1;
			right.pressed = true;
			return true;
		} else if (evt.key.keysym.sym == SDLK_w) {
			up.downs += 1;
			up.pressed = true;
			return true;
		} else if (evt.key.keysym.sym == SDLK_s) {
			down.downs += 1;
			down.pressed = true;
			return true;
		}
		else if (evt.key.keysym.sym == SDLK_r) {
			reset();
			return true;
		}
		else if (evt.key.keysym.sym == SDLK_1) {
			reset(1);
			return true;
		}
		else if (evt.key.keysym.sym == SDLK_2) {
			reset(2);
			return true;
		}
		else if (evt.key.keysym.sym == SDLK_3) {
			reset(3);
			return true;
		}
		else if (evt.key.keysym.sym == SDLK_4) {
			reset(4);
			return true;
		}
		else if (evt.key.keysym.sym == SDLK_5) {
			reset(5);
			return true;
		}
		else if (evt.key.keysym.sym == SDLK_6) {
			reset(6);
			return true;
		}
		else if (evt.key.keysym.sym == SDLK_7) {
			reset(7);
			return true;
		}
		else if (evt.key.keysym.sym == SDLK_8) {
			reset(8);
			return true;
		}
		else if (evt.key.keysym.sym == SDLK_9) {
			reset(9);
			return true;
		}
	} else if (evt.type == SDL_KEYUP) {
		if (evt.key.keysym.sym == SDLK_a) {
			left.pressed = false;
			return true;
		} else if (evt.key.keysym.sym == SDLK_d) {
			right.pressed = false;
			return true;
		} else if (evt.key.keysym.sym == SDLK_w) {
			up.pressed = false;
			return true;
		} else if (evt.key.keysym.sym == SDLK_s) {
			down.pressed = false;
			return true;
		}
	}

	return false;
}

void PlayMode::update(float elapsed) {
	// Move target based on player input
	glm::vec3 d_target(0, 0, 0);
	if (left.pressed) {
		d_target[0] --;
	}
	if (right.pressed) {
		d_target[0] ++;
	}
	if (up.pressed) {
		d_target[1] ++;
	}
	if (down.pressed) {
		d_target[1] --;
	}
	if (glm::length(d_target) > 1) {
		d_target /= glm::length(d_target);
	}
	target += d_target * target_move_speed * elapsed;

	if (glm::length(target) > max_target_distance) {
		target *= max_target_distance / glm::length(target);
	}
	
	// Update rotation of laser cannon
	cannon_theta = atan2(target.y, target.x) - (float)M_PI / 2.f;
	glm::vec3 cannon_pos = cannon->make_local_to_world() * glm::vec4(0.f, 0.f, 0.f, 1.f);
	cannon_phi = atan2(-cannon_pos[2], glm::length(target));

	// Rotate the cannon
	cannon->rotation = glm::angleAxis(cannon_phi, glm::vec3(1.f, 0.f, 0.f));
	pole->rotation = glm::angleAxis(cannon_theta, glm::vec3(0.f, 0.f, 1.f));

	// Update the rotation of the laser beam
	laser_offset = (laser_offset + laser_spin_speed * elapsed);
	laser_offset -= (int)(laser_offset / (2.f * M_PI)) * 2.f * (float)M_PI;
	

	// Destroy the tile targeted by the laser
	GroundTile* targeted_tile = getTileAtCoords(glm::vec2(target[0], target[1]));
	if (targeted_tile != nullptr) {
		destroyTile(targeted_tile);
	}

	// Update tiles
	for (size_t i = 0; i < map.size(); i ++) {
		GroundTile* tile = &map[i];
		if (isInRing(tile) && tile->exists) {
			if (tile->falling) {
				tile->velocity += gravity * elapsed;
				tile->transform->position += tile->velocity * elapsed;
				if (tile->transform->position[2] < -100) {
					tile->exists = false;
				}
			}
		}
	}

	// Update meeples
	bool all_sick = true;
	GroundTile* root = nullptr;
	for (size_t i = 0; i < meeples.size(); i++) {
		Meeple* meeple = meeples[i];
		
		if (meeple->tile->exists) {
			if (meeple->sick) {
				// If this meeple is sick and we don't already have one, set root tile to this one
				if (root == nullptr) {
					root = meeple->tile;
				}
			}
			else {
				// If this meeple isn't sick, they aren't all sick
				all_sick = false;
			}

			// Set jump constants based on meeple's health
			float wait = jump_wait_duration;
			float h = jump_height;
			float dur = jump_duration;
			if (meeple->sick) {
				wait *= sick_wait_multiplier;
				h *= sick_height_multiplier;
				dur *= sick_duration_multiplier;
			}

			// Update jumps from tile to tile
			if (!meeple->jumping_to) {
				// If the meeple is not jumping, it sits on its tile
				meeple->drawable->transform->position = meeple->tile->transform->position;

				// Update timer until next jump
				meeple->jump_wait_timer += elapsed;
				
				if (meeple->tile->exists && !meeple->tile->falling && meeple->jump_wait_timer >= wait) {
					// Reset wait timer
					meeple->jump_wait_timer -= wait;

					// Get vector of neighbors
					glm::ivec2 pos = getTilePos(meeple->tile);
					std::vector<GroundTile*> neighbors;
					std::vector<float> weights;
					glm::vec2 bounce = glm::vec2(0.f, 0.f);
					for (size_t j = 0; j < directions.size(); j++) {
						GroundTile* neighbor = getTileAtPos(pos + directions[j]);
						if (neighbor == nullptr) {
							bounce -= (glm::vec2)directions[j];
						}
						else if (neighbor->meeple == nullptr) {
							neighbors.push_back(neighbor);
							weights.push_back(powf(1.f + glm::dot(meeple->dir, (glm::vec2)directions[j]), 2));
						}
					}

					// Jump to random neighbor
					if (neighbors.size() > 0) {
						// Normalize weights
						float total_weight = 0;
						for (size_t j = 0; j < weights.size(); j++) {
							total_weight += weights[j];
						}
						for (size_t j = 0; j < weights.size(); j++) {
							weights[j] /= total_weight;
						}

						// Get weighted random neighbor index to jump to
						float r = (float)rand() / (float)RAND_MAX;
						float cumulative_weight = 0;
						size_t index = 0;
						for (size_t j = 0; j < weights.size(); j++) {
							cumulative_weight += weights[j];
							if (cumulative_weight >= r) {
								index = j;
								glm::vec2 new_dir = getTileCoords(neighbors[j]) - getTileCoords(meeple->tile);
								meeple->dir += new_dir * inertia;
								meeple->dir += bounce * elasticity;
								if (glm::length(meeple->dir) > 0) {
									meeple->dir /= glm::length(meeple->dir);
								}
								else {
									meeple->dir = glm::vec2(0.f, 1.f);
								}
								break;
							}
						}

						// Update state of meeple and relevant tiles
						meeple->jump_timer = 0;
						meeple->jumping_to = neighbors[index];
						meeple->jumping_from = meeple->tile;
						meeple->tile->meeple = nullptr;
						meeple->tile = meeple->jumping_to;
						meeple->tile->meeple = meeple;
					}
				}
			}
			else {
				// If the meeple is jumping, update its position along a parabolic arc
				meeple->jump_timer += elapsed;
				if (meeple->jump_timer <= dur) {
					glm::vec2 p0 = getTileCoords(meeple->jumping_from);
					glm::vec2 p1 = getTileCoords(meeple->jumping_to);
					glm::vec2 delta = p1 - p0;
					glm::vec2 p = p0 + (meeple->jump_timer / dur) * delta;

					// Based on z = a (t + b)^2 + c, where t is normalized between 0 and 1
					float a = -4 * h;
					float b = -0.5f;
					float c = h;
					float z = a * powf(meeple->jump_timer / dur + b, 2) + c;
					meeple->drawable->transform->position = glm::vec3(p[0], p[1], z);
				}
				else {
					// End jump
					meeple->drawable->transform->position = meeple->tile->transform->position;
					meeple->jumping_to = nullptr;
					meeple->jumping_from = nullptr;
				}
			}

			// If sick, update infection counters of nearby meeples
			if (meeple->sick && !meeple->tile->falling) {
				glm::ivec2 pos = getTilePos(meeple->tile);
				for (int y = -infection_radius; y <= infection_radius; y++) {
					for (int x = -infection_radius; x <= infection_radius; x++) {
						GroundTile* tile = getTileAtPos(pos + glm::ivec2(x, y));
						if (tile && tile->meeple) {
							infect(tile->meeple, elapsed);
						}
					}
				}
			}
		}
		else {
			// Hide meeple if its tile has been destroyed
			meeple->drawable->transform->position[2] = -100.f;
		}
	}

	// If all meeples are sick, the game is over
	if (all_sick && root != nullptr) {
		dropChunk(root);
	}

	// Create cloud
	cloud_timer -= elapsed;
	if (cloud_timer <= 0) {
		Cloud* cloud = new Cloud;
		cloud->transform = new Scene::Transform();
		scene.drawables.emplace_back(cloud->transform);
		cloud->drawable = &scene.drawables.back();
		setMesh(cloud->drawable, "Cloud");

		auto randFloat = [](float min, float max) -> float {
			return min + (float)rand() / (float)RAND_MAX * (max - min);
		};

		// Set cloud position
		cloud->transform->position.x = randFloat(-map_radius * 1.5, map_radius * 1.5);
		cloud->transform->position.y = camera->transform->position.y;
		cloud->transform->position.z = randFloat(cloud_min_altitude, cloud_max_altitude);
		
		// Set cloud size
		cloud->transform->scale.x = randFloat(cloud_min_width, cloud_max_width);
		cloud->transform->scale.y = randFloat(cloud_min_length, cloud_max_length);
		cloud->transform->scale.z = randFloat(cloud_min_height, cloud_max_height);

		// Set cloud speed
		cloud->speed = randFloat(cloud_min_speed, cloud_max_speed);

		clouds.push_back(cloud);

		cloud_timer = randFloat(cloud_min_interval, cloud_max_interval);
	}

	// Update clouds
	for (size_t i = 0; i < clouds.size(); i++) {
		clouds[i]->transform->position += glm::vec3(0.f, clouds[i]->speed * elapsed, 0.f);
		if (clouds[i]->transform->position.y - clouds[i]->transform->scale.y > 100.f) {
			deleteCloud(i);
			i--;
		}
	}

	//reset button press counters:
	left.downs = 0;
	right.downs = 0;
	up.downs = 0;
	down.downs = 0;
}

void PlayMode::draw(glm::uvec2 const &drawable_size) {
	//update camera aspect ratio for drawable:
	camera->aspect = float(drawable_size.x) / float(drawable_size.y);

	//set up light type and position for lit_color_texture_program:
	glUseProgram(lit_color_texture_program->program);
	glUniform1i(lit_color_texture_program->LIGHT_TYPE_int, 1);
	glUniform3fv(lit_color_texture_program->LIGHT_DIRECTION_vec3, 1, glm::value_ptr(glm::vec3(0.f, 1.f,-.5f)));
	glUniform3fv(lit_color_texture_program->LIGHT_ENERGY_vec3, 1, glm::value_ptr(glm::vec3(1.0f, 1.0f, 0.95f)));
	glUseProgram(0);

	glClearColor(0.f, 0.5f, 1.f, 1.0f);
	glClearDepth(1.0f); //1.0 is actually the default value to clear the depth buffer to, but FYI you can change it.
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS); //this is the default depth comparison function, but FYI you can change it.

	GL_ERRORS(); //print any errors produced by this setup code

	scene.draw(*camera);

	{ // Draw laser beam
		glm::mat4 world_to_clip = camera->make_projection() * glm::mat4(camera->transform->make_world_to_local());
		DrawLines lines(world_to_clip);

		// Convert hsv triple to u8vec4 rgb(a), with a always 255
		auto hsvToRgb = [](float h, uint8_t s, uint8_t v) -> glm::u8vec4 {
			h -= (int)floor(h / 360) * 360;
			glm::tvec3<double> rgb = glm::rgbColor(glm::tvec3<double>((double)h, (double)s / 100, (double)v / 100));
			return glm::u8vec4(rgb[0] * 255, rgb[1] * 255, rgb[2] * 255, 255);
		};

		// Draw beam
		for (float theta = 0 + laser_offset; theta < 2.f * M_PI + laser_offset; theta += 2.f * (float)M_PI / num_beams) {
			glm::vec3 a = nozzle->make_local_to_world() * glm::vec4(0.f + cosf(theta) * laser_radius, 0.f, 0.f + sinf(theta) * laser_radius, 1.f);
			glm::vec3 b = target + nozzle->make_local_to_world() * glm::vec4(cosf(theta) * laser_radius / 2, 0.f, sinf(theta) * laser_radius / 2, 1.f) - (nozzle->make_local_to_world() * glm::vec4(0.f, 0.f, 0.f, 1.f));
			lines.draw(a, b, hsvToRgb(laser_hue_amplitude * sin((theta - laser_offset) * laser_hue_frequency) + laser_hue_offset, 100, 100));
		}

		// Draw sparks where the laser hits the ground plane
		for (size_t i = 0; i < num_sparks; i++) {
			glm::vec3 a = target;
			auto rng = []() -> float {
				return 2 * (float)rand() / (float)RAND_MAX - 1;
			};
			glm::vec3 b = a + glm::vec3(rng() * spark_size, rng() * spark_size, rng() * spark_size);
			lines.draw(a, b, hsvToRgb(laser_hue_amplitude * rng() + laser_hue_offset, 100, 100));
		}
	}

	{ // Print text
		glDisable(GL_DEPTH_TEST);

		float aspect = float(drawable_size.x) / float(drawable_size.y);
		DrawLines lines(glm::mat4(
			1.0f / aspect, 0.0f, 0.0f, 0.0f,
			0.0f, 1.0f, 0.0f, 0.0f,
			0.0f, 0.0f, 1.0f, 0.0f,
			0.0f, 0.0f, 0.0f, 1.0f
		));

		float H = 0.1f;
		float ofs = 2.0f / drawable_size.y;
		auto drawText = [&](std::string str, float x, float y, glm::u8vec4 color, bool bold = false) {
			int thickness = 1;
			if (bold) {
				thickness = 3;
			}
			for (int i = -thickness; i <= thickness; i++) {
				for (int j = -thickness; j <= thickness; j++) {
					lines.draw_text(str,
						glm::vec3(-aspect + x * H + ofs * j, -1.0 + y * H + ofs * i, 0.0),
						glm::vec3(H, 0.0f, 0.0f), glm::vec3(0.0f, H, 0.0f),
						glm::u8vec4(0x00, 0x00, 0x00, 0x00));
				}
			}

			if (bold) {
				for (int i = -1; i <= 1; i++) {
					for (int j = -1; j <= 1; j++) {
						lines.draw_text(str,
							glm::vec3(-aspect + x * H + ofs * j, -1.0 + y * H + ofs * i, 0.0),
							glm::vec3(H, 0.0f, 0.0f), glm::vec3(0.0f, H, 0.0f),
							color);
					}
				}
			}
			else {
				lines.draw_text(str,
					glm::vec3(-aspect + x * H, -1.0 + y * H, 0.0),
					glm::vec3(H, 0.0f, 0.0f), glm::vec3(0.0f, H, 0.0f),
					color);
			}
		};
		
		// Print score text
		drawText("Score: ", 0.2f, 2.f, glm::u8vec4(0xff, 0xff, 0xff, 0xff));
		glm::u8vec4 score_color = score >= 0 ? glm::u8vec4(0xff, 0xff, 0xff, 0xff) : glm::u8vec4(0xd8, 0x00, 0x00, 0xff);
		drawText(std::to_string(score), 3.f, 2.f, score_color, true);
		
		// Print grade text
		drawText("Grade: ", 0.2f, 0.2f, glm::u8vec4(0xff, 0xff, 0xff, 0xff));
		std::string grade;
		glm::u8vec4 grade_color;

		// Calculate the grade and associated color
		float grade_score;
		if (perfect_score == 0) {
			// This happens if the player hasn't done anything yet
			if (score == 0) {
				grade = "";
				grade_color = glm::u8vec4(0x00, 0x00, 0x00, 0xff);
			}
			else {
				// Now that I think about it, I'm pretty sure this should never happen
				grade = "F";
				grade_color = glm::u8vec4(0x40, 0x00, 0x00, 0xff);
			}
		}
		else {
			// Grade is score as percentage of perfect score
			grade_score = (float)score / (float)perfect_score;

			// Calculate grade and color in normal case
			// Note that we kind of fudge the grades here; S is 90+, A is 80+, etc. It just feels better.
			if (grade_score >= 0.9f) {
				grade_color = glm::u8vec4(0x00, 0xd8, 0xff, 0xff);
				if (grade_score > 0.99999f) {
					grade = "SS";
				}
				else if (grade_score >= 0.97f) {
					grade = "S+";
				}
				else if (grade_score >= 0.93f) {
					grade = "S";
				}
				else {
					grade = "S-";
				}
			}
			else if (grade_score >= 0.8f) {
				grade_color = glm::u8vec4(0x00, 0xd8, 0x00, 0xff);
				if (grade_score >= 0.87f) {
					grade = "A+";
				}
				else if (grade_score >= 0.83f) {
					grade = "A";
				}
				else {
					grade = "A-";
				}
			}
			else if (grade_score >= 0.7f) {
				grade_color = glm::u8vec4(0xff, 0xd8, 0x00, 0xff);
				if (grade_score >= 0.77f) {
					grade = "B+";
				}
				else if (grade_score >= 0.73f) {
					grade = "B";
				}
				else {
					grade = "B-";
				}
			}
			else if (grade_score >= 0.6f) {
				grade_color = glm::u8vec4(0xff, 0x80, 0x00, 0xff);
				if (grade_score >= 0.67f) {
					grade = "C+";
				}
				else if (grade_score >= 0.63) {
					grade = "C";
				}
				else {
					grade = "C-";
				}
			}
			else if (grade_score >= 0.5f) {
				grade_color = glm::u8vec4(0xd8, 0x00, 0x00, 0xff);
				if (grade_score >= 0.57f) {
					grade = "D+";
				}
				else if (grade_score >= 0.53) {
					grade = "D";
				}
				else {
					grade = "D-";
				}
			}
			else {
				grade = "F";
				grade_color = glm::u8vec4(0x40, 0x00, 0x00, 0xff);
			}
		}
		
		// Print the grade
		drawText(grade, 3.f, .2f, grade_color, true);

		// Print some brief instructions
		drawText("R: Restart. 1-9: Change Difficulty.", 14.f, .2f, glm::u8vec4(0xFF, 0xFF, 0xFF, 0xFF));
	}
}
