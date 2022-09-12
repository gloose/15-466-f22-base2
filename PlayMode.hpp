#include "Mode.hpp"

#include "Scene.hpp"

#include <glm/glm.hpp>

#include <vector>
#include <deque>
#include <array>

struct PlayMode : Mode {
	PlayMode();
	virtual ~PlayMode();

	//functions called by main loop:
	virtual bool handle_event(SDL_Event const &, glm::uvec2 const &window_size) override;
	virtual void update(float elapsed) override;
	virtual void draw(glm::uvec2 const &drawable_size) override;

	//----- game state -----

	//input tracking:
	struct Button {
		uint8_t downs = 0;
		uint8_t pressed = 0;
	} left, right, down, up;

	//local copy of the game scene (so code can change it during gameplay):
	Scene scene;

	// Camera
	Scene::Camera* camera = nullptr;

	// Important parts of the cannon
	Scene::Transform* pole = nullptr;
	Scene::Transform* cannon = nullptr;
	Scene::Transform* nozzle = nullptr;

	// Describes the spinning of the laser beam (purely aesthetic)
	float laser_offset = 0.f;
	float laser_spin_speed = 8.f;
	float laser_radius = 0.1f;
	float laser_hue_amplitude = 40.f;
	float laser_hue_frequency = 3.f;
	float laser_hue_offset = 0;
	uint8_t num_beams = 20;

	// Describes rotation angle of the cannon
	float cannon_theta = 0.f;
	float cannon_phi = 0.f;

	// Pre-define meeple struct to allow for circular dependency
	struct Meeple;

	// Struct representing a single cube of the ground plane
	struct GroundTile {
		Scene::Transform* transform = nullptr;
		bool falling = false;
		glm::vec3 velocity = glm::vec3(0.f, 0.f, 0.f);
		bool exists = false;
		bool visited = false;
		GroundTile* visited_for = nullptr;
		uint8_t chunk = 0;
		Meeple* meeple = nullptr;
	};

	// Struct representing a meeple
	struct Meeple {
		Scene::Drawable* drawable = nullptr;
		GroundTile* tile = nullptr;
		GroundTile* jumping_to = nullptr;
		GroundTile* jumping_from = nullptr;
		float jump_timer = 0;
		bool sick = false;
		float infection = 0;
		glm::vec2 dir = glm::vec2(0.f, 0.f);
		float jump_wait_timer = 0;
	};

	// Descibes meeple jumps
	const float jump_duration = 0.5f;
	const float jump_height = 0.5f;
	const float jump_wait_duration = 0.5f;
	const float sick_duration_multiplier = .25f;
	const float sick_height_multiplier = 1.f;
	const float sick_wait_multiplier = .25f;

	// Map stored as a grid of tiles
	static const int map_radius = 24;
	static const int hole_radius = 4;
	static const int map_diameter = 2 * map_radius + 1;
	const float tile_size = .5f;
	std::array<GroundTile, map_diameter * map_diameter> map;

	// Vector of all meeples on the map
	std::vector<Meeple*> meeples;

	// Acceleration due to gravity
	glm::vec3 gravity = glm::vec3(0.f, 0.f, -9.8f);

	// Describes the target point
	glm::vec3 target = glm::vec3(0.f, -1.5f, 0.f);
	size_t num_sparks = 128;
	float spark_size = .2f;

	// Arrays of direction vectors, used for checking tile neighbors iteratively
	std::array<glm::ivec2, 4> directions = {
		glm::ivec2(1, 0),
		glm::ivec2(-1, 0),
		glm::ivec2(0, 1),
		glm::ivec2(0, -1)
	};

	// Describes spread of infection between meeples
	const float infection_threshold = 60.f;
	const size_t num_thresholds = 4;
	const int infection_radius = 4;

	// Describes the degree to which meeples maintain their direction, or "bounce" off of a gap they can't cross
	const float inertia = 0.5f;
	const float elasticity = 1.f;

	// Fraction of tiles holding a meeple at the start of the game
	const float population_density = 0.05f;

	// Score earned this game
	int score = 0;

	// Score values for various actions
	int score_healthy = 100;
	int score_sick = 10;
	int score_kill = -200;

	// Speed at which the ground target moves
	const float target_move_speed = 4.f;

	// Helper functions
	bool isInRing(glm::ivec2 pos);
	bool isInRing(GroundTile* tile);
	glm::ivec2 getTilePos(GroundTile* tile);
	glm::vec2 getTileCoords(GroundTile* tile);
	GroundTile* getTileAtPos(glm::ivec2 pos);
	GroundTile* getTileAtCoords(glm::vec2 coords);
	void destroyTile(GroundTile* tile);
	void setMesh(Scene::Drawable* drawable, std::string mesh);
	void infect(Meeple* meeple, float exposure);
	void reset();
};
